from pathlib import Path
from ase.io import read, write
from ase import Atom
from ase.data import covalent_radii
from ase.constraints import FixAtoms
from ase.neighborlist import neighbor_list
import numpy as np
import json

def _surface_oxygen_indices(slab, surface_threshold: float):
    positions = slab.get_positions()
    symbols = slab.get_chemical_symbols()
    z_top = positions[:, 2].max()
    surface_O = [
        i for i, s in enumerate(symbols)
        if s == "O" and (z_top - positions[i, 2]) < surface_threshold
    ]
    return surface_O, float(z_top)

def _nearest_k_oxygen_sites(slab, center_o_index: int, candidates: list[int], k: int = 8):
    """Return list of (o_index, distance_A) sorted by distance to center oxygen."""
    center_pos = slab.positions[center_o_index]
    dists = []
    for idx in candidates:
        if idx == center_o_index:
            continue
        if slab[idx].symbol != "O":
            continue
        dist = float(np.linalg.norm(slab.positions[idx] - center_pos))
        dists.append((int(idx), dist))
    dists.sort(key=lambda x: x[1])
    return dists[:k]


def _one_hot_bits(n: int, hot_index: int) -> str:
    bits = ["0"] * n
    bits[hot_index] = "1"
    return "".join(bits)

def setup_structure(
    input_file: str = "inputs/slab_clean_2x2_canonical.in",
    o_index: int = 20,
    oh_distance: float = 1.0,
    z_freeze: float = 21.01,
    neighbor_cutoff: float = 1.5,
    surface_threshold: float = 1.5,  # New param for surface detection
    outputs_dir="outputs",
):
    # Create output directory
    outdir = Path(outputs_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Adjust H covalent radius
    covalent_radii[1] = 0.6
    
    print(f"Reading {input_file}...")
    slab = read(input_file, format="espresso-in")
    
    # Disable PBC for xTB
    slab.set_pbc((False, False, False))
    
    print(f"Initial structure: {len(slab)} atoms")
    print(f"Cell: {slab.cell.cellpar()}")
    
    # Find surface oxygen atoms
    positions = slab.get_positions()
    symbols = slab.get_chemical_symbols()
    z_top = positions[:, 2].max()
    
    surface_O = [
        i for i, s in enumerate(symbols)
        if s == "O" and (z_top - positions[i, 2]) < surface_threshold  # Use param
    ]
    
    print(f"\nSurface O atoms (within {surface_threshold} Å of top): {surface_O}")
    
    if o_index not in surface_O:
        print(f"  Warning: O[{o_index}] is not in surface list!")
        print(f"   O[{o_index}] is at z={positions[o_index, 2]:.2f} Å")
        print(f"   Top of slab: z={z_top:.2f} Å")
    
    # Add hydrogen atom
    O_pos = slab.positions[o_index]
    H_pos = O_pos + np.array([0.0, 0.0, oh_distance])
    slab.append(Atom("H", position=H_pos))
    h_index = len(slab) - 1
    
    print(f"\nAdded H atom:")
    print(f"  Index: {h_index}")
    print(f"  Position: {H_pos}")
    print(f"  Total atoms: {len(slab)}")
    
    # Freeze bottom layers
    z = slab.positions[:, 2]
    freeze_mask = z < z_freeze
    slab.set_constraint(FixAtoms(mask=freeze_mask))
    n_frozen = int(freeze_mask.sum())
    n_mobile = len(slab) - n_frozen
    
    print(f"\nFreezing atoms below z={z_freeze:.1f} Å:")
    print(f"  Frozen: {n_frozen} atoms")
    print(f"  Mobile: {n_mobile} atoms")
    
    # Check neighbors around H
    i, j, d = neighbor_list("ijd", slab, cutoff=neighbor_cutoff)  # Simplified cutoff
    neighbors = [(int(jj), float(dd)) for ii, jj, dd in zip(i, j, d) if ii == h_index]
    
    print(f"\nNeighbors within {neighbor_cutoff} Å of H:")
    for jj, dd in neighbors:
        sym = slab[jj].symbol
        pos = slab[jj].position
        print(f"  {sym}[{jj}]: {dd:.3f} Å at ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})")

    # Build metadata dict
    meta = {
        "input_file": str(input_file),
        "o_index": int(o_index),
        "h_index": int(h_index),
        "oh_distance_A": float(oh_distance),
        "o_position_A": [float(x) for x in O_pos],
        "h_position_A": [float(x) for x in H_pos],
        "zmax_A": float(z_top),  # Fixed: was zmax (undefined)
        "z_freeze_A": float(z_freeze),
        "neighbor_cutoff_A": float(neighbor_cutoff),
        "neighbors_of_H": [{"index": int(jj), "distance_A": float(dd)} for jj, dd in neighbors],
        "num_atoms_total": int(len(slab)),
    }

    # Save prepared structure
    (outdir / "structures").mkdir(parents=True, exist_ok=True)
    (outdir / "results").mkdir(parents=True, exist_ok=True)

    output_file = outdir / "structures" / f"slab_H_o{o_index}_ready.traj"
    write(str(output_file), slab)

    meta_file = outdir / "results" / f"metadata_o{o_index}.json"
    meta_file.write_text(json.dumps(meta, indent=2))
    
    print(f"\n{'='*60}")
    print(f" Structure prepared successfully!")
    print(f"{'='*60}")
    print(f"Output files:")
    print(f"  Structure: {output_file}")
    print(f"  Metadata:  {meta_file}")
    print(f"\nNext step:")
    print(f"  python run_optimization.py {output_file}")
    print(f"{'='*60}")
    
    return slab, meta

def setup_batch_nearby_oxygen_sites(
    input_file: str = "inputs/slab_clean_2x2_canonical.in",
    center_o_index: int = 20,
    k: int = 8,
    oh_distance: float = 1.0,
    z_freeze: float = 21.01,
    neighbor_cutoff: float = 1.5,
    surface_threshold: float = 1.5,
    outputs_dir: str = "outputs",
    include_center: bool = False,
):
    """
    Batch-generate H* structures by placing H on the top-k nearest SURFACE oxygen sites
    around center_o_index. Writes each case into outputs_dir/batch_centeroXX/oYYY/{structures,results}.
    """
    outdir = Path(outputs_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read slab once to find nearest sites (no H added here)
    slab0 = read(input_file, format="espresso-in")
    slab0.set_pbc((False, False, False))

    surface_O, z_top = _surface_oxygen_indices(slab0, surface_threshold=surface_threshold)
    if center_o_index not in surface_O:
        print(f"Warning: center O[{center_o_index}] is not a surface O by threshold {surface_threshold} Å.")
        print(f"Surface O list: {surface_O}")

    # Candidate target O sites: surface O atoms (optionally include the center)
    candidates = list(surface_O)
    if not include_center and center_o_index in candidates:
        candidates.remove(center_o_index)

    nearest = _nearest_k_oxygen_sites(slab0, center_o_index=center_o_index, candidates=candidates, k=k)

    batch_root = outdir / f"batch_centerO{center_o_index}"
    batch_root.mkdir(parents=True, exist_ok=True)

    print(f"\nNearest {k} surface O sites to center O[{center_o_index}] (by Cartesian distance):")
    for rank, (o_idx, dist) in enumerate(nearest):
        print(f"  rank {rank}: O[{o_idx}] at {dist:.3f} Å")

    # Build each structure by calling existing setup_structure with per-site outputs
    for rank, (target_o_idx, dist) in enumerate(nearest):
        run_dir = batch_root / f"o{target_o_idx}"
        run_dir.mkdir(parents=True, exist_ok=True)

        slab, meta = setup_structure(
            input_file=input_file,
            o_index=int(target_o_idx),
            oh_distance=oh_distance,
            z_freeze=z_freeze,
            neighbor_cutoff=neighbor_cutoff,
            surface_threshold=surface_threshold,
            outputs_dir=str(run_dir),
        )

        # Add batch metadata + 1-hot bits label for Hamiltonian fitting later
        meta["center_o_index"] = int(center_o_index)
        meta["target_o_index"] = int(target_o_idx)
        meta["site_rank"] = int(rank)
        meta["distance_to_center_O_A"] = float(dist)
        meta["k_sites"] = int(k)
        meta["bits"] = _one_hot_bits(k, rank)  # single-occupancy (1-hot) encoding

        # Overwrite metadata file with augmented fields
        meta_file = Path(run_dir) / "results" / f"metadata_o{target_o_idx}.json"
        meta_file.write_text(json.dumps(meta, indent=2))

    print(f"\nBatch complete. Outputs in: {batch_root}")
    print(f"Next step: run your optimization script for each prepared .traj in each folder.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Prepare H* adsorption structure(s)")
    parser.add_argument("--input", default="inputs/slab_clean_2x2_canonical.in", help="Input structure file")

    # Single-run args (existing)
    parser.add_argument("--o-index", type=int, default=20, help="Oxygen atom index for H placement")
    parser.add_argument("--oh-dist", type=float, default=1.0, help="O-H distance (Å)")
    parser.add_argument("--z-freeze", type=float, default=20.0, help="Freeze atoms below this z (Å)")
    parser.add_argument("--neighbor-cutoff", type=float, default=1.5, help="Neighbor cutoff (Å)")
    parser.add_argument("--surface-threshold", type=float, default=1.5, help="Threshold for surface O detection (Å)")
    parser.add_argument("--outputs", default="outputs", help="Outputs directory")

    # Batch args (new)
    parser.add_argument("--batch-nearby-k", type=int, default=0,
                        help="If >0, prepare a batch over the k nearest surface O sites around --center-o-index.")
    parser.add_argument("--center-o-index", type=int, default=20,
                        help="Center oxygen index used to rank nearby O sites (batch mode).")
    parser.add_argument("--include-center", action="store_true",
                        help="Include the center oxygen itself as a candidate target site (batch mode).")

    args = parser.parse_args()

    if args.batch_nearby_k and args.batch_nearby_k > 0:
        setup_batch_nearby_oxygen_sites(
            input_file=args.input,
            center_o_index=args.center_o_index,
            k=args.batch_nearby_k,
            oh_distance=args.oh_dist,
            z_freeze=args.z_freeze,
            neighbor_cutoff=args.neighbor_cutoff,
            surface_threshold=args.surface_threshold,
            outputs_dir=args.outputs,
            include_center=args.include_center,
        )
    else:
        slab, meta = setup_structure(
            input_file=args.input,
            o_index=args.o_index,
            oh_distance=args.oh_dist,
            z_freeze=args.z_freeze,
            neighbor_cutoff=args.neighbor_cutoff,
            surface_threshold=args.surface_threshold,
            outputs_dir=args.outputs
        )

