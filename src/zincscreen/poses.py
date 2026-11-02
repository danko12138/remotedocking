"""Retain every finite energy; report geometric flags separately from scoring."""
from collections import Counter
from pathlib import Path
import math
import re

from .chemistry import pdbqt_atoms


def vina_poses(path, expected_atoms):
    text = Path(path).read_text()
    blocks = re.findall(r"^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", text, re.M | re.S)
    if not blocks or [int(x[0]) for x in blocks] != list(range(1, len(blocks) + 1)):
        raise ValueError("Empty/truncated/nonsequential Vina pose output")
    expected = {r["serial"]: r["type"] for r in expected_atoms}
    poses = []
    for number, block in blocks:
        energy = re.findall(r"^REMARK VINA RESULT:\s*(\S+)", block, re.M)
        if len(energy) != 1 or not math.isfinite(float(energy[0])):
            raise ValueError("Missing/nonfinite pose energy")
        atoms = pdbqt_atoms(block, text=True)
        if {r["serial"]: r["type"] for r in atoms} != expected:
            raise ValueError("Pose atom identity/count changed")
        poses.append({"pose": int(number), "energy_kcal_mol": float(energy[0]),
                      "xyz_by_serial": {r["serial"]: r["xyz"] for r in atoms}})
    return poses


def gnina_poses(path, ligand_sdf):
    from rdkit import Chem
    source = Chem.RemoveHs(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)[0])
    poses = []
    for i, mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=False), 1):
        if mol is None or not mol.HasProp("minimizedAffinity"):
            raise ValueError("Missing/invalid GNINA SDF record or minimizedAffinity")
        energy = float(mol.GetProp("minimizedAffinity"))
        if not math.isfinite(energy):
            raise ValueError("Nonfinite GNINA energy")
        heavy = Chem.RemoveHs(mol)
        order = heavy.GetSubstructMatch(source, useChirality=True)
        if heavy.GetNumAtoms() != source.GetNumAtoms() or len(order) != source.GetNumAtoms():
            raise ValueError("GNINA output ligand graph/stereochemistry mismatch")
        xyz = heavy.GetConformer().GetPositions()
        if not all(math.isfinite(float(v)) for row in xyz for v in row):
            raise ValueError("Nonfinite GNINA coordinates")
        row = {"pose": i, "energy_kcal_mol": energy, "source_heavy_xyz": xyz[list(order)].tolist()}
        for key in ("CNNscore", "CNNaffinity"):
            if mol.HasProp(key):
                value = float(mol.GetProp(key))
                if not math.isfinite(value):
                    raise ValueError("Nonfinite " + key)
                row[key] = value
        poses.append(row)
    if not poses:
        raise ValueError("Empty GNINA output")
    return poses


def geometry_flags(pose, ligand_sdf, mapping, receptor, box):
    from rdkit import Chem
    import numpy as np
    mol = Chem.RemoveHs(Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)[0])
    if "source_heavy_xyz" in pose:
        xyz = np.asarray(pose.pop("source_heavy_xyz"))
    else:
        coords = pose.pop("xyz_by_serial")
        xyz = np.array([coords[s] for s in mapping["source_heavy_to_pdbqt_serial"]])
    pt = Chem.GetPeriodicTable()
    invalid_bonds = 0
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        radii = pt.GetRcovalent(a.GetAtomicNum()) + pt.GetRcovalent(b.GetAtomicNum())
        length = math.dist(xyz[a.GetIdx()], xyz[b.GetIdx()])
        invalid_bonds += not .55 * radii <= length <= 1.6 * radii
    rec = pdbqt_atoms(receptor)
    rxyz = np.array([r["xyz"] for r in rec])
    distances = np.linalg.norm(xyz[:, None] - rxyz[None, :], axis=2)
    minimum = float(distances.min())
    center, size = np.array(box["center"]), np.array(box["size"])
    outside = int(np.any(np.abs(xyz - center) > size / 2 + 1.0, axis=1).sum())
    # These are proximity contacts, not inferred hydrogen bonds or binding evidence.
    contact_serials = {rec[i]["serial"] for i in np.flatnonzero(np.any(distances <= 4.0, axis=0))}
    residues = set()
    for line in Path(receptor).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in contact_serials:
            residues.add(line[17:27].strip())
    return {"invalid_bonds": invalid_bonds, "minimum_receptor_distance_A": minimum,
            "heavy_atoms_outside_box_plus_1A": outside,
            "qc_eligible": invalid_bonds == 0 and minimum >= 1.2,
            "contact_residues_within_4A": sorted(residues),
            "qc_limitation": "Coarse geometry/contact screen; visual review and detailed chemistry checks remain required"}
