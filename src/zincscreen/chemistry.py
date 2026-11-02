"""Explicit input state, bounded conformer generation, strict native conversion."""
from pathlib import Path
import math

from .common import sha, save


def geometry_ok(mol):
    from rdkit import Chem
    if mol is None or mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        return False
    xyz = mol.GetConformer().GetPositions()
    if not all(math.isfinite(float(v)) for row in xyz for v in row):
        return False
    pt = Chem.GetPeriodicTable()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        radii = pt.GetRcovalent(a.GetAtomicNum()) + pt.GetRcovalent(b.GetAtomicNum())
        if not .55 * radii <= math.dist(xyz[a.GetIdx()], xyz[b.GetIdx()]) <= 1.6 * radii:
            return False
    return True


def stereo_check(source, mol):
    from rdkit import Chem
    expected, actual = Chem.Mol(source), Chem.Mol(mol)
    Chem.SanitizeMol(expected)
    Chem.SanitizeMol(actual)
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    Chem.AssignStereochemistryFrom3D(actual, replaceExistingTags=True)
    Chem.AssignStereochemistry(actual, cleanIt=True, force=True)
    for atom in expected.GetAtoms():
        if atom.HasProp("_CIPCode"):
            other = actual.GetAtomWithIdx(atom.GetIdx())
            if not other.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != other.GetProp("_CIPCode"):
                raise ValueError("Specified stereochemistry changed")
    for bond in expected.GetBonds():
        if str(bond.GetStereo()) in ("STEREOE", "STEREOZ", "STEREOCIS", "STEREOTRANS"):
            other = actual.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if str(other.GetStereo()) != str(bond.GetStereo()):
                raise ValueError("Specified double-bond stereochemistry changed")


def pdbqt_atoms(path_or_text, *, text=False):
    content = path_or_text if text else Path(path_or_text).read_text()
    rows = []
    for line in content.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            typ = line.split()[-1]
            xyz = [float(line[i:i + 8]) for i in (30, 38, 46)]
            charge = float(line[70:76])
            if not all(math.isfinite(v) for v in (*xyz, charge)):
                raise ValueError("Nonfinite PDBQT coordinates or charges")
            if typ not in ("H", "HD", "HS"):
                rows.append({"serial": int(line[6:11]), "type": typ, "xyz": xyz})
    if not rows or len({r["serial"] for r in rows}) != len(rows):
        raise ValueError("Empty or duplicated heavy-atom serials")
    return rows


def native_mapping(sdf, pdbqt):
    from rdkit import Chem
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if not geometry_ok(mol):
        raise ValueError("Unusable SDF geometry")
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    atoms = pdbqt_atoms(pdbqt)
    if len(atoms) != len(heavy):
        raise ValueError("Native heavy-atom count mismatch")
    elements = {"A": "C", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "SA": "S"}
    xyz = mol.GetConformer().GetPositions()[[a.GetIdx() for a in heavy]]
    q = np.array([r["xyz"] for r in atoms])
    distances = np.linalg.norm(xyz[:, None] - q[None, :], axis=2)
    mismatch = np.array([[a.GetSymbol() != elements.get(r["type"], r["type"]) for r in atoms] for a in heavy])
    ii, jj = linear_sum_assignment(distances + mismatch * 1e6)
    if mismatch[ii, jj].any() or distances[ii, jj].max(initial=0) > .002:
        raise ValueError("Native conversion changed elements or geometry")
    serials = [None] * len(heavy)
    for i, j in zip(ii, jj):
        serials[int(i)] = atoms[int(j)]["serial"]
    return {"source_heavy_to_pdbqt_serial": serials,
            "maximum_displacement_A": float(distances[ii, jj].max(initial=0)),
            "sdf_sha256": sha(sdf), "pdbqt_sha256": sha(pdbqt)}


def prepare_one(record, folder, seed=42):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    from openbabel import openbabel as ob
    source = Chem.MolFromSmiles(record["smiles"])
    if source is None:
        raise ValueError("SMILES_PARSE_FAILED")
    if len(Chem.GetMolFrags(source)) != 1 or any(a.GetAtomicNum() == 0 for a in source.GetAtoms()):
        raise ValueError("MULTICOMPONENT_OR_NONCONCRETE_INPUT: no automatic salt stripping")
    key = Chem.MolToInchiKey(source)
    if not key or record.get("inchikey", key) != key:
        raise ValueError("SOURCE_IDENTITY_MISMATCH")
    attempts, mol = [], None
    for chosen in (seed, seed + 1, seed + 2):
        trial = Chem.AddHs(source)
        params = AllChem.ETKDGv3()
        params.randomSeed, params.numThreads, params.maxIterations = chosen, 1, 1000
        params.enforceChirality = True
        params.useRandomCoords = chosen != seed
        code = AllChem.EmbedMolecule(trial, params)
        if code < 0:
            attempts.append({"seed": chosen, "status": "EMBED_FAILED"})
            continue
        field, convergence = None, None
        if AllChem.MMFFHasAllMoleculeParams(trial):
            field = "MMFF94s"
            convergence = AllChem.MMFFOptimizeMolecule(trial, mmffVariant=field, maxIters=1000)
        elif AllChem.UFFHasAllMoleculeParams(trial):
            field = "UFF"
            convergence = AllChem.UFFOptimizeMolecule(trial, maxIters=1000)
        attempts.append({"seed": chosen, "field": field, "optimization_code": convergence})
        # New library inputs require a supported, converged force field.
        if field is None or convergence != 0 or not geometry_ok(trial):
            continue
        Chem.SanitizeMol(trial)
        try:
            stereo_check(source, trial)
        except ValueError as exc:
            attempts[-1]["error"] = str(exc)
            continue
        mol = trial
        break
    if mol is None:
        raise ValueError("NO_VALID_CONVERGED_CONFORMER: " + str(attempts))
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    sdf, pdbqt = folder / "ligand.sdf", folder / "ligand.pdbqt"
    mol.SetProp("_Name", record["id"])
    mol.SetProp("SOURCE_SMILES", record["smiles"])
    mol.SetProp("SOURCE_INCHIKEY", key)
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(mol)
    reader = ob.OBConversion()
    reader.SetInFormat("sdf")
    native = ob.OBMol()
    if not reader.ReadFile(native, str(sdf)):
        raise ValueError("SDF_READ_FAILED")
    # Strict RDKit parameter handling avoids silently zeroing unknown elements.
    # Transfer only after confirming unchanged SDF atom order/elements/coordinates.
    AllChem.ComputeGasteigerCharges(mol, nIter=12, throwOnParamFailure=True)
    charges = [float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()]
    charge = Chem.GetFormalCharge(source)
    if not all(math.isfinite(v) for v in charges) or abs(sum(charges) - charge) > .05:
        raise ValueError("PARTIAL_CHARGE_SUM_MISMATCH")
    native_atoms = list(ob.OBMolAtomIter(native))
    if len(native_atoms) != mol.GetNumAtoms():
        raise ValueError("SDF_CONVERSION_ATOM_COUNT_CHANGED")
    xyz = mol.GetConformer().GetPositions()
    for i, (atom, q) in enumerate(zip(native_atoms, charges)):
        if atom.GetAtomicNum() != mol.GetAtomWithIdx(i).GetAtomicNum() or math.dist(xyz[i], (atom.GetX(), atom.GetY(), atom.GetZ())) > .002:
            raise ValueError("SDF_CONVERSION_ATOM_ORDER_CHANGED")
        if atom.GetFormalCharge() != mol.GetAtomWithIdx(i).GetFormalCharge():
            raise ValueError("SDF_CONVERSION_FORMAL_CHARGE_CHANGED")
        atom.SetPartialCharge(q)
    native.SetPartialChargesPerceived(True)
    writer = ob.OBConversion()
    writer.SetOutFormat("pdbqt")
    writer.AddOption("p", ob.OBConversion.OUTOPTIONS)  # Preserve input atom indices.
    if not writer.WriteFile(native, str(pdbqt)):
        raise ValueError("PDBQT_WRITE_FAILED")
    writer.CloseOutFile()
    rounded = sum(float(line[70:76]) for line in pdbqt.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM")))
    if abs(rounded - charge) > .1:
        raise ValueError("PDBQT_CHARGE_SUM_MISMATCH")
    mapping = native_mapping(sdf, pdbqt)
    save(folder / "mapping.json", mapping)
    return {"status": "READY", "inchikey": key, "formal_charge": charge,
            "method": "ETKDGv3+converged-forcefield+RDKit-strict-Gasteiger+OpenBabel-PDBQT", "attempts": attempts,
            "charge_method": {"name": "RDKit Gasteiger", "nIter": 12, "throwOnParamFailure": True,
                              "partial_charge_sum": sum(charges), "pdbqt_charge_sum": rounded},
            "rdkit_version": rdBase.rdkitVersion,
            "unspecified_source_stereo": any(str(s.specified) == "Unspecified" for s in Chem.FindPotentialStereo(source)),
            "limitation": "One representative supplied chemical state; no pH/tautomer ensemble. Unspecified stereo is not experimentally assigned.",
            "files": {name: sha(folder / name) for name in ("ligand.sdf", "ligand.pdbqt", "mapping.json")}}
