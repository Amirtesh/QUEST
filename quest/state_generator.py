"""
QUEST – Module 2: The State Generator
=======================================
Enumerates biologically relevant tautomers of a ligand, builds a 3-D
conformer for each, and runs a CREST conformer/minima search to obtain
the QM-level minimum-energy geometry for every protonation/tautomer state.

Workflow
--------
1. generate_tautomers  – SMILES → ranked RDKit Mol list via
                          TautomerEnumerator.
2. embed_3d            – add H atoms and generate an ETKDGv3 conformer.
3. run_crest_minimum_search – subprocess wrapper around the CREST binary,
                              returning contents of ``crest_best.xyz``.
4. update_mol_from_xyz_string – transplant QM coordinates back onto the
                                RDKit Mol object.
5. prep_ligand_states  – orchestrator: tautomers → embed → CREST → SDF.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Geometry import rdGeometry

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Function 1
# ---------------------------------------------------------------------------

def generate_tautomers(
    smiles: str,
    max_tautomers: int = 3,
) -> list[Chem.Mol]:
    """Enumerate the top tautomers of a molecule from its SMILES.

    Uses RDKit's ``TautomerEnumerator`` which ranks tautomers by a scoring
    function that penalises unstable forms (e.g. enols, vinyl amines).

    Parameters
    ----------
    smiles:
        Input SMILES string.
    max_tautomers:
        Maximum number of tautomers to return (highest-scoring first).

    Returns
    -------
    list[Chem.Mol]
        Up to *max_tautomers* sanitized RDKit Mol objects.  The canonical
        tautomer (index 0) is always included.

    Raises
    ------
    ValueError
        If the SMILES cannot be parsed or sanitization fails.
    """
    mol: Chem.Mol | None = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            f"RDKit could not parse SMILES: '{smiles}'. "
            "Check for typos or non-standard notation."
        )

    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException as exc:
        raise ValueError(
            f"Sanitization failed for SMILES '{smiles}': {exc}"
        ) from exc

    enumerator = rdMolStandardize.TautomerEnumerator()
    tautomers: list[Chem.Mol] = list(enumerator.Enumerate(mol))

    if not tautomers:
        logger.warning(
            "TautomerEnumerator returned 0 tautomers for '%s'; using input mol.",
            smiles,
        )
        tautomers = [mol]

    selected = tautomers[:max_tautomers]
    logger.info(
        "Generated %d tautomer(s) for '%s' (requested max %d).",
        len(selected),
        smiles,
        max_tautomers,
    )
    return selected


# ---------------------------------------------------------------------------
# Function 2
# ---------------------------------------------------------------------------

def embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """Add explicit hydrogens and generate a 3-D conformer via ETKDGv3.

    Parameters
    ----------
    mol:
        An RDKit Mol object (no conformer required).

    Returns
    -------
    Chem.Mol
        The same molecule with H atoms added and a single embedded 3-D
        conformer.

    Raises
    ------
    RuntimeError
        If ``EmbedMolecule`` returns ``-1`` (embedding failed because the
        distance geometry could not be satisfied – rare but possible for
        strained or exotic scaffolds).
    """
    mol_h: Chem.Mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42          # reproducibility
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True

    result: int = AllChem.EmbedMolecule(mol_h, params)
    if result == -1:
        smiles_repr = Chem.MolToSmiles(mol)
        raise RuntimeError(
            f"ETKDGv3 embedding failed for molecule: {smiles_repr}. "
            "The distance-geometry solver could not place all atoms. "
            "Consider simplifying the structure or pre-generating a 2-D layout."
        )

    logger.info(
        "3-D conformer generated successfully (%d atoms).", mol_h.GetNumAtoms()
    )
    return mol_h


# ---------------------------------------------------------------------------
# Function 3
# ---------------------------------------------------------------------------

def run_crest_minimum_search(
    xyz_path: str,
    threads: int = 4,
) -> str:
    """Run a CREST iMTD-GC conformer search and return the best geometry.

    The search is run inside a temporary scratch directory to avoid littering
    the working directory with CREST output files.  After completion the
    contents of ``crest_best.xyz`` are read and returned as a string.

    Parameters
    ----------
    xyz_path:
        Absolute (or relative) path to the input XYZ file.
    threads:
        Number of CPU threads to pass to CREST via ``-T``.

    Returns
    -------
    str
        Full contents of ``crest_best.xyz`` – the lowest-energy geometry
        found by CREST, ready to be parsed by
        :func:`update_mol_from_xyz_string`.

    Raises
    ------
    RuntimeError
        If CREST exits with a non-zero return code, exceeds the 2-hour
        timeout, or does not produce ``crest_best.xyz``.
    """
    xyz_abs: str = os.path.abspath(xyz_path)

    cmd: list[str] = [
        "crest", xyz_abs,
        "--gfn2",
        "--alpb", "water",
        "--quick",
        "--ewin", "3.0",
        "-T", str(threads),
    ]

    logger.info("Launching CREST: %s", " ".join(cmd))

    # Prevent OpenBLAS from spawning its own OpenMP threads, which conflict
    # with CREST's internal threading model and cause hangs on AMD/Ryzen.
    # CREST manages parallelism itself via -T; OpenBLAS must stay serial.
    crest_env = os.environ.copy()
    crest_env["OMP_NUM_THREADS"] = "1"
    crest_env["MKL_NUM_THREADS"] = "1"
    crest_env["OPENBLAS_NUM_THREADS"] = "1"

    with tempfile.TemporaryDirectory() as scratch_dir:
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=scratch_dir,
                env=crest_env,
            )

            import threading, sys

            stdout_lines: list[str] = []

            def _stream() -> None:
                for line in proc.stdout:
                    stdout_lines.append(line)
                    print(line, end="", flush=True)

            reader = threading.Thread(target=_stream, daemon=True)
            reader.start()

            try:
                proc.wait(timeout=7200)
            except subprocess.TimeoutExpired:
                proc.kill()
                reader.join(timeout=5)
                raise RuntimeError(
                    f"CREST timed out after 7200 s for input: {xyz_path}. "
                    "Consider reducing the molecule size or increasing thread count."
                )

            reader.join(timeout=30)
            stdout_text = "".join(stdout_lines)

            if proc.returncode != 0:
                raise RuntimeError(
                    f"CREST returned non-zero exit code {proc.returncode} "
                    f"for input: {xyz_path}.\n"
                    f"--- CREST output (last 40 lines) ---\n"
                    + "\n".join(stdout_text.splitlines()[-40:])
                )

            logger.info(
                "CREST finished successfully for %s", os.path.basename(xyz_abs)
            )

        except RuntimeError:
            raise
        except Exception as exc:
            raise RuntimeError(
                f"Unexpected error running CREST for {xyz_path}: {exc}"
            ) from exc

        # --- read crest_best.xyz before the tempdir is deleted ---------------
        best_xyz_path = os.path.join(scratch_dir, "crest_best.xyz")
        if not os.path.isfile(best_xyz_path):
            raise RuntimeError(
                f"CREST completed but 'crest_best.xyz' was not found in {scratch_dir}. "
                f"CREST stdout tail:\n" + "\n".join(stdout_text.splitlines()[-20:])
            )

        with open(best_xyz_path, "r", encoding="utf-8") as fh:
            best_xyz_contents: str = fh.read()

    logger.info("crest_best.xyz read successfully (%d chars).", len(best_xyz_contents))
    return best_xyz_contents


# ---------------------------------------------------------------------------
# Function 4
# ---------------------------------------------------------------------------

def update_mol_from_xyz_string(mol: Chem.Mol, xyz_string: str) -> Chem.Mol:
    """Overwrite the conformer coordinates of an RDKit Mol from an XYZ string.

    CREST preserves the original atom ordering produced by the input XYZ
    file, so a direct positional mapping (atom index ↔ XYZ row) is safe.

    XYZ format assumed::

        <n_atoms>
        <comment line>
        El  x  y  z
        El  x  y  z
        ...

    Parameters
    ----------
    mol:
        RDKit Mol with an existing conformer (e.g. from :func:`embed_3d`).
    xyz_string:
        String contents of a valid XYZ file whose atom count matches
        ``mol.GetNumAtoms()``.

    Returns
    -------
    Chem.Mol
        The mol object with conformer positions updated to the CREST
        QM-optimised coordinates.

    Raises
    ------
    ValueError
        If the atom count in the XYZ string does not match the molecule,
        or if the XYZ block cannot be parsed.
    """
    lines = xyz_string.strip().splitlines()

    # --- parse header --------------------------------------------------------
    try:
        n_atoms_xyz: int = int(lines[0].strip())
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"Could not parse atom count from XYZ header: {lines[:2]}"
        ) from exc

    n_atoms_mol: int = mol.GetNumAtoms()
    if n_atoms_xyz != n_atoms_mol:
        raise ValueError(
            f"Atom count mismatch: XYZ has {n_atoms_xyz} atoms but molecule "
            f"has {n_atoms_mol} atoms. Ensure CREST input was generated from "
            f"the same molecule."
        )

    # --- parse coordinate block (lines[2:] skips n_atoms + comment) ----------
    coord_lines = lines[2:]
    if len(coord_lines) < n_atoms_mol:
        raise ValueError(
            f"XYZ string has only {len(coord_lines)} coordinate lines "
            f"but {n_atoms_mol} are required."
        )

    conf: Chem.Conformer = mol.GetConformer()
    for i, line in enumerate(coord_lines[:n_atoms_mol]):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(
                f"XYZ line {i + 3} is malformed (expected 'El x y z'): '{line}'"
            )
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError as exc:
            raise ValueError(
                f"Could not parse coordinates on XYZ line {i + 3}: '{line}'"
            ) from exc

        conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))

    logger.info(
        "Updated %d atom positions from CREST XYZ output.", n_atoms_mol
    )
    return mol


# ---------------------------------------------------------------------------
# Function 5
# ---------------------------------------------------------------------------

def prep_ligand_states(
    smiles: str,
    output_dir: str = ".",
    max_tautomers: int = 3,
    threads: int = 4,
) -> list[str]:
    """Enumerate tautomers, embed each in 3-D, run CREST, and save as SDF.

    This is the top-level orchestrator for Module 2.  For each tautomer the
    pipeline is:

    .. code-block:: text

        SMILES
          └─ generate_tautomers()
               └─ [for each tautomer i]
                    embed_3d()
                    → write temp .xyz
                    → run_crest_minimum_search()
                    → update_mol_from_xyz_string()
                    → save tautomer_{i}_qm.sdf

    Parameters
    ----------
    smiles:
        Input SMILES string of the ligand.
    output_dir:
        Directory where per-tautomer SDF files are written.
    max_tautomers:
        Maximum number of tautomers to enumerate (passed to
        :func:`generate_tautomers`).
    threads:
        CPU threads to pass to CREST via ``-T``.

    Returns
    -------
    list[str]
        Absolute paths to the saved SDF files, one per tautomer processed
        successfully.

    Raises
    ------
    ValueError
        If the SMILES is invalid.
    RuntimeError
        If embedding or CREST fails for all tautomers.
    """
    os.makedirs(output_dir, exist_ok=True)
    tautomers: list[Chem.Mol] = generate_tautomers(smiles, max_tautomers=max_tautomers)

    saved_paths: list[str] = []

    for i, taut_mol in enumerate(tautomers):
        taut_smiles = Chem.MolToSmiles(taut_mol)
        logger.info("Processing tautomer %d: %s", i, taut_smiles)

        # Step a – embed 3-D conformer
        try:
            mol_3d: Chem.Mol = embed_3d(taut_mol)
        except RuntimeError as exc:
            logger.warning(
                "Skipping tautomer %d – 3-D embedding failed: %s", i, exc
            )
            continue

        # Step b – write to a temporary XYZ file
        with tempfile.NamedTemporaryFile(
            suffix=".xyz",
            delete=False,
            mode="w",
            encoding="utf-8",
        ) as tmp_xyz:
            Chem.MolToXYZFile(mol_3d, tmp_xyz.name)
            tmp_xyz_path: str = tmp_xyz.name

        try:
            # Step c – CREST minimum search
            logger.info(
                "Running CREST on tautomer %d (input: %s) …", i, tmp_xyz_path
            )
            best_xyz_string: str = run_crest_minimum_search(
                tmp_xyz_path, threads=threads
            )

            # Step d – update mol with QM-optimised coordinates
            mol_qm: Chem.Mol = update_mol_from_xyz_string(mol_3d, best_xyz_string)

        except RuntimeError as exc:
            logger.warning(
                "Skipping tautomer %d – CREST failed: %s", i, exc
            )
            continue
        except ValueError as exc:
            logger.warning(
                "Skipping tautomer %d – coordinate update failed: %s", i, exc
            )
            continue
        finally:
            # Clean up the temporary xyz input file
            try:
                os.unlink(tmp_xyz_path)
            except OSError:
                pass

        # Step e – save as SDF
        sdf_name = f"tautomer_{i}_qm.sdf"
        sdf_path = os.path.abspath(os.path.join(output_dir, sdf_name))
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol_qm)
        writer.close()

        saved_paths.append(sdf_path)
        logger.info("Saved %s", sdf_path)

    if not saved_paths:
        raise RuntimeError(
            "All tautomers failed processing. "
            "Check SMILES validity, embedding logs, and CREST installation."
        )

    logger.info(
        "prep_ligand_states complete: %d SDF file(s) written.", len(saved_paths)
    )
    return saved_paths
