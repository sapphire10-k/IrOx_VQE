#!/bin/bash
set -euo pipefail

SITE="69"

# --- Inputs
INPUT_SLAB="inputs/slab_clean_2x2_canonical.in"

# --- Where we put the prepared H* state (match your workflow convention)
STATE_DIR="outputs/states/H_star"
READY_TRAJ="${STATE_DIR}/structures/slab_H_o${SITE}_ready.traj"

# --- Where optimization writes (your convention differs across scripts; we keep it simple)
OPT_OUTDIR="${STATE_DIR}/results"
# Try to use a predictable name for the optimized traj; if your optimizer uses a different one,
# we’ll discover it after the run.
OPT_TRAJ_EXPECTED="${OPT_OUTDIR}/slab_H_o${SITE}_ready_opt.traj"

# --- Downstream outputs
CLUSTER_OUTDIR="outputs/clusters"
ORBITAL_OUTDIR="outputs/orbital_analysis/neutral"
HAM_OUTDIR="outputs/hamiltonians"

# --- Parameters (choose conservative defaults)
OH_DIST="1.0"
Z_FREEZE="21.01"
SURFACE_THRESHOLD="1.5"
NEIGHBOR_CUTOFF="1.5"

CLUSTER_CUTOFF="6.0"

# Orbital analysis canonical case (neutral, closed-shell)
CHARGE="0"
UHF="0"
SHELL_CUTOFF="0.10"
FRONTIER_WINDOW_EV="3.0"
TOPN="12"

# Hamiltonian (keep light for CI; tune as you like)
HAM_CHARGE="0"
HAM_SPIN="0"
BASIS="sto-3g"
METHOD="RHF"
REGION_CUTOFF="6.0"
N_OCC="6"
N_VIRT="6"

# --- Runner helper: prefer micromamba if present, else conda
run_env () {
  local ENV_NAME="$1"
  shift
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n "$ENV_NAME" "$@"
  elif command -v conda >/dev/null 2>&1; then
    conda run -n "$ENV_NAME" "$@"
  else
    echo "ERROR: Neither micromamba nor conda found on PATH."
    exit 1
  fi
}

echo "==> Canonical example: slab setup -> optimize -> cluster -> orbitals -> Hamiltonian (o${SITE})"

mkdir -p "${STATE_DIR}/structures" "${STATE_DIR}/results" "${CLUSTER_OUTDIR}" "${ORBITAL_OUTDIR}" "${HAM_OUTDIR}"

echo "==> Step 1: Slab setup (place H on O index ${SITE})"
# IMPORTANT: this script is in your repo as scripts/iro2_slab_setup.py (not the uploaded filename)
run_env iro2 python -u scripts/iro2_slab_setup.py \
  --input "${INPUT_SLAB}" \
  --o-index "${SITE}" \
  --oh-dist "${OH_DIST}" \
  --z-freeze "${Z_FREEZE}" \
  --neighbor-cutoff "${NEIGHBOR_CUTOFF}" \
  --surface-threshold "${SURFACE_THRESHOLD}" \
  --outputs "${STATE_DIR}"

test -f "${READY_TRAJ}" || (echo "ERROR: setup did not produce ${READY_TRAJ}" && exit 1)
echo "Prepared: ${READY_TRAJ}"

echo "==> Step 2: Geometry optimisation"
# Your CI uses: optimization_iro2_test.py --input ... --output ...
run_env iro2 python -u scripts/optimization_iro2_test.py \
  "${READY_TRAJ}" \
  --output "${OPT_OUTDIR}"

# Locate optimised traj
OPT_TRAJ=""
if [[ -f "${OPT_TRAJ_EXPECTED}" ]]; then
  OPT_TRAJ="${OPT_TRAJ_EXPECTED}"
else
  OPT_TRAJ="$(ls -1 "${OPT_OUTDIR}"/*.traj 2>/dev/null | head -n 1 || true)"
fi

if [[ -z "${OPT_TRAJ}" || ! -f "${OPT_TRAJ}" ]]; then
  echo "ERROR: No optimised .traj found in ${OPT_OUTDIR}"
  echo "Contents:"
  ls -la "${OPT_OUTDIR}" || true
  exit 1
fi
echo "Optimised structure: ${OPT_TRAJ}"

echo "==> Step 3: Extract cluster"
run_env iro2 python -u scripts/extract_cluster.py \
  --traj "${OPT_TRAJ}" \
  --center_atom_index "${SITE}" \
  --cutoff "${CLUSTER_CUTOFF}" \
  --outdir "${CLUSTER_OUTDIR}" \
  --tag "o${SITE}"

echo "Clusters directory:"
ls -la "${CLUSTER_OUTDIR}" || true

# Locate the cluster structure file produced (robust search)
CLUSTER_PATH="$(find "${CLUSTER_OUTDIR}" -type f \
  \( -name "*o${SITE}*.traj" -o -name "*o${SITE}*.xyz" -o -name "*o${SITE}*.pdb" \) \
  | head -n 1 || true)"

if [[ -z "${CLUSTER_PATH}" || ! -f "${CLUSTER_PATH}" ]]; then
  echo "ERROR: Could not locate cluster file for o${SITE} in ${CLUSTER_OUTDIR}"
  echo "Debug (first 50 files):"
  find "${CLUSTER_OUTDIR}" -type f | head -n 50 || true
  exit 1
fi
echo "Using cluster: ${CLUSTER_PATH}"

echo "==> Step 4: Orbital analysis (xTB cluster)"
run_env iro2 python scripts/orbital_analysis_xtb_cluster.py \
  "${CLUSTER_PATH}" \
  --site "${SITE}" \
  --charge "${CHARGE}" \
  --uhf "${UHF}" \
  --shell-cutoff "${SHELL_CUTOFF}" \
  --frontier-window-ev "${FRONTIER_WINDOW_EV}" \
  --topn "${TOPN}" \
  --outdir "${ORBITAL_OUTDIR}"

echo "Orbital analysis outputs:"
ls -la "${ORBITAL_OUTDIR}" || true

#echo "==> Step 5: Build fermionic Hamiltonian (PySCF) in ham env"
#run_env ham python scripts/build_fermionic_hamiltonian_pyscf.py \
  # "${CLUSTER_PATH}" \
  # --site "${SITE}" \
  # --charge "${HAM_CHARGE}" \
  # --spin "${HAM_SPIN}" \
  # --basis "${BASIS}" \
  # --method "${METHOD}" \
  # --region-cutoff "${REGION_CUTOFF}" \
  # --n-occ "${N_OCC}" \
  # --n-virt "${N_VIRT}" \
  # --outdir "${HAM_OUTDIR}"

echo "==> Step 5: Build fermionic Hamiltonian (optional, heavy-element)"

echo "Skipping Hamiltonian build in canonical example."
echo "Reason: heavy-element (Ir) basis/ECP setup is required and not suitable for lightweight CI."
echo "See scripts/build_fermionic_hamiltonian_pyscf.py for full functionality."

#echo "Hamiltonian outputs:"
#ls -la "${HAM_OUTDIR}" || true

echo "==> Example workflow complete."
