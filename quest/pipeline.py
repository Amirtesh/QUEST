"""
QUEST – Master Orchestrator: The Pipeline
==========================================
Combines all three modules into a single end-to-end run:

  Module 2 → tautomer/QM states
  Module 3 → docking (P2Rank + Vina)
  Module 1 → conformational strain per docked pose

Output is a tidy Pandas DataFrame (and CSV) with one row per pose, ranked
by Vina affinity and annotated with a Viable flag based on strain cutoff.

Usage (programmatic)
--------------------
    from quest.pipeline import run_quest_pipeline
    df = run_quest_pipeline(
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        receptor_pdb="receptor.pdb",
        receptor_pdbqt="receptor.pdbqt",
        size=(20.0, 20.0, 20.0),
        output_dir="./quest_run",
    )
    print(df)

Usage (CLI)
-----------
    quest run "CC(=O)Oc1ccccc1C(=O)O" receptor.pdb receptor.pdbqt \\
        --size 20 20 20 --output-dir ./quest_run -v
"""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem

from quest.state_generator import prep_ligand_states
from quest.docking_engine import dock_ligand
from quest.strain_evaluator import calculate_strain

logger = logging.getLogger(__name__)

# Strain cutoff (kcal mol⁻¹) above which a pose is flagged as non-viable
STRAIN_CUTOFF: float = 7.0


def _extract_vina_affinity(pose: Chem.Mol) -> float | None:
    """Extract the Vina binding affinity from an RDKit Mol produced by OpenBabel.

    OpenBabel encodes the Vina score in one of several properties depending on
    the obabel version.  We try them in priority order:

    1. ``minimizedAffinity``  – set by recent obabel (most reliable)
    2. ``REMARK``             – raw REMARK block; parse first float after
                               "VINA RESULT:" or the first numeric token
    3. ``_Name`` / title     – some obabel versions embed the score in the
                               molecule title

    Parameters
    ----------
    pose:
        An RDKit Mol object read from a Vina-output SDF.

    Returns
    -------
    float | None
        Binding affinity in kcal mol⁻¹, or ``None`` if it cannot be parsed.
    """
    # 1 – minimizedAffinity (most common)
    for prop in ("minimizedAffinity", "minimizedaffinity", "Affinity", "affinity"):
        if pose.HasProp(prop):
            try:
                return float(pose.GetProp(prop))
            except (ValueError, TypeError):
                pass

    # 2 – REMARK block
    if pose.HasProp("REMARK"):
        remark = pose.GetProp("REMARK")
        for line in remark.splitlines():
            line_upper = line.upper()
            if "VINA RESULT" in line_upper or "AFFINITY" in line_upper:
                tokens = line.split()
                for token in tokens:
                    try:
                        return float(token)
                    except ValueError:
                        continue
        # no labelled line – try any float in the whole REMARK
        for token in remark.split():
            try:
                return float(token)
            except ValueError:
                continue

    # 3 – molecule title / _Name
    if pose.HasProp("_Name"):
        title = pose.GetProp("_Name")
        for token in title.split():
            try:
                return float(token)
            except ValueError:
                continue

    return None


def run_quest_pipeline(
    smiles: str,
    receptor_pdb: str,
    receptor_pdbqt: str,
    size: tuple[float, float, float],
    output_csv: str = "quest_results.csv",
    output_dir: str = ".",
    strain_cutoff: float = STRAIN_CUTOFF,
    max_tautomers: int = 3,
    threads: int = 4,
    exhaustiveness: int = 16,
) -> pd.DataFrame:
    """Run the full QUEST pipeline and return a results DataFrame.

    Steps
    -----
    1. Enumerate tautomers and generate QM-minimum geometries (Module 2).
    2. Dock each tautomer into the receptor (Module 3).
    3. For every docked pose, compute the conformational strain (Module 1).
    4. Assemble results, add ``Viable`` flag, save CSV, return DataFrame.

    Parameters
    ----------
    smiles:
        Input SMILES of the ligand.
    receptor_pdb:
        Receptor PDB file for P2Rank pocket detection.
    receptor_pdbqt:
        Receptor PDBQT file for AutoDock Vina.
    size:
        ``(sx, sy, sz)`` – Vina docking box dimensions in Å.
    output_csv:
        Path where the results CSV is saved.
    output_dir:
        Root directory for all intermediate and final output files.
    strain_cutoff:
        Strain threshold (kcal mol⁻¹) for the ``Viable`` flag
        (default: 7.0 kcal mol⁻¹).
    max_tautomers:
        Maximum tautomers to generate and process.
    threads:
        CPU threads for CREST (``-T``).
    exhaustiveness:
        AutoDock Vina exhaustiveness.

    Returns
    -------
    pd.DataFrame
        One row per docked pose with columns:
        ``Tautomer``, ``Pose``, ``Vina_Affinity``, ``QM_Strain``, ``Viable``.
        Sorted by ``Vina_Affinity`` (ascending, most negative first), then
        ``QM_Strain`` (ascending).

    Raises
    ------
    RuntimeError
        If Module 2 fails entirely (no tautomers could be processed).
    ValueError
        If the SMILES is invalid.
    """
    os.makedirs(output_dir, exist_ok=True)
    results: list[dict] = []

    # ------------------------------------------------------------------
    # Stage 1: Tautomer enumeration + QM geometry optimisation (Module 2)
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("STAGE 1 – State generation for: %s", smiles)
    logger.info("=" * 60)

    tautomer_sdfs: list[str] = prep_ligand_states(
        smiles,
        output_dir=output_dir,
        max_tautomers=max_tautomers,
        threads=threads,
    )
    logger.info("State generator produced %d SDF file(s).", len(tautomer_sdfs))

    # ------------------------------------------------------------------
    # Stage 2 & 3: Dock each tautomer, then score each pose
    # ------------------------------------------------------------------
    for taut_idx, tautomer_qm_sdf in enumerate(tautomer_sdfs):
        taut_label = os.path.basename(tautomer_qm_sdf)

        logger.info("-" * 60)
        logger.info(
            "STAGE 2 – Docking tautomer %d/%d: %s",
            taut_idx + 1, len(tautomer_sdfs), taut_label,
        )
        logger.info("-" * 60)

        # --- dock --------------------------------------------------------
        try:
            docked_sdf_path: str = dock_ligand(
                qm_sdf=tautomer_qm_sdf,
                receptor_pdbqt=receptor_pdbqt,
                size=size,
                p2rank_pdb=receptor_pdb,
                output_dir=output_dir,
                exhaustiveness=exhaustiveness,
            )
        except (RuntimeError, FileNotFoundError, ValueError) as exc:
            logger.warning(
                "Docking failed for tautomer %d (%s) – skipping. Error: %s",
                taut_idx, taut_label, exc,
            )
            continue

        logger.info("Docked SDF: %s", docked_sdf_path)

        # --- iterate over poses ------------------------------------------
        supplier = Chem.SDMolSupplier(docked_sdf_path, sanitize=False)
        n_poses = len(supplier)
        logger.info("Reading %d pose(s) from docked SDF.", n_poses)

        for pose_idx, pose in enumerate(supplier):
            if pose is None:
                logger.warning(
                    "Pose %d of %s could not be read by RDKit – skipping.",
                    pose_idx, taut_label,
                )
                continue

            # --- Vina affinity -------------------------------------------
            affinity: float | None = _extract_vina_affinity(pose)
            if affinity is None:
                logger.warning(
                    "Could not extract Vina affinity for tautomer %d pose %d.",
                    taut_idx, pose_idx,
                )

            # --- strain: write single pose to temp SDF and call Module 1 --
            strain: float | None = None
            tmp_pose_sdf: str | None = None
            try:
                with tempfile.NamedTemporaryFile(
                    suffix=".sdf",
                    delete=False,
                    dir=output_dir,
                    mode="w",
                    encoding="utf-8",
                ) as tmp_f:
                    tmp_pose_sdf = tmp_f.name

                writer = Chem.SDWriter(tmp_pose_sdf)
                writer.write(pose)
                writer.close()

                logger.info(
                    "STAGE 3 – Strain calc: tautomer %d pose %d …",
                    taut_idx, pose_idx,
                )
                strain = calculate_strain(
                    ref_sdf=tautomer_qm_sdf,
                    docked_sdf=tmp_pose_sdf,
                )
                logger.info(
                    "  Vina affinity = %s kcal mol⁻¹ | Strain = %.4f kcal mol⁻¹",
                    f"{affinity:.3f}" if affinity is not None else "N/A",
                    strain,
                )

            except (RuntimeError, ValueError) as exc:
                logger.warning(
                    "Strain calculation failed for tautomer %d pose %d: %s",
                    taut_idx, pose_idx, exc,
                )
                strain = None
            finally:
                if tmp_pose_sdf and os.path.isfile(tmp_pose_sdf):
                    try:
                        os.unlink(tmp_pose_sdf)
                    except OSError:
                        pass

            results.append({
                "Tautomer": taut_label,
                "Pose": pose_idx,
                "Vina_Affinity": affinity,
                "QM_Strain": strain,
            })

    # ------------------------------------------------------------------
    # Stage 4: Build DataFrame, add Viable flag, save CSV
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("STAGE 4 – Assembling results (%d pose entries).", len(results))
    logger.info("=" * 60)

    df = pd.DataFrame(results, columns=["Tautomer", "Pose", "Vina_Affinity", "QM_Strain"])

    # Viable: True if strain is known AND below the cutoff
    df["Viable"] = df["QM_Strain"].apply(
        lambda s: bool(s is not None and s <= strain_cutoff)
        if s is not None else False
    )

    # Sort: best Vina score first (most negative), then lowest strain
    df.sort_values(
        by=["Vina_Affinity", "QM_Strain"],
        ascending=[True, True],
        na_position="last",
        inplace=True,
    )
    df.reset_index(drop=True, inplace=True)

    output_csv_abs = os.path.abspath(output_csv)
    df.to_csv(output_csv_abs, index=False)
    logger.info("Results saved → %s", output_csv_abs)

    return df


def run_quest_batch(
    input_path: str,
    receptor_pdb: str,
    receptor_pdbqt: str,
    size: tuple[float, float, float],
    output_csv: str = "quest_batch_results.csv",
    output_dir: str = ".",
    **kwargs,
) -> pd.DataFrame:
    """Run :func:`run_quest_pipeline` over a collection of ligands.

    Accepts either a text file of SMILES strings or a directory of SDF files
    and processes each ligand independently, collecting all pose results into a
    single master DataFrame.

    Parameters
    ----------
    input_path:
        One of:

        * A ``.txt`` or ``.smi`` file where each non-empty, non-comment line
          contains a SMILES string optionally followed by a whitespace-separated
          molecule ID, e.g. ``CC(=O)O aspirin``.
        * A directory path.  Every ``.sdf`` file found directly inside the
          directory is read; the stem of each filename is used as the
          molecule ID.
    receptor_pdb:
        Receptor PDB file for P2Rank pocket detection.
    receptor_pdbqt:
        Receptor PDBQT file for AutoDock Vina.
    size:
        ``(sx, sy, sz)`` – Vina docking box dimensions in Å.
    output_csv:
        Path where the master results CSV (all ligands combined) is saved.
    output_dir:
        Root directory under which a per-molecule sub-directory is created
        for each ligand's intermediate and final output files.
    **kwargs:
        Any additional keyword arguments forwarded verbatim to
        :func:`run_quest_pipeline` (e.g. ``strain_cutoff``, ``max_tautomers``,
        ``threads``, ``exhaustiveness``).

    Returns
    -------
    pd.DataFrame
        Concatenation of all per-ligand DataFrames, with an extra
        ``Ligand_ID`` column prepended.  Returns an empty DataFrame if every
        ligand failed.

    Raises
    ------
    ValueError
        If *input_path* is neither a recognised file type nor a directory.
    """
    input_p = Path(input_path)

    # ------------------------------------------------------------------
    # 1. Parse input into a list of (smiles, mol_id) tuples
    # ------------------------------------------------------------------
    molecules: list[tuple[str, str]] = []

    if input_p.is_file() and input_p.suffix.lower() in (".txt", ".smi"):
        logger.info("Batch input: SMILES file %s", input_p)
        with input_p.open("r", encoding="utf-8") as fh:
            for lineno, raw_line in enumerate(fh, start=1):
                line = raw_line.strip()
                # skip blank lines and comment lines
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                smiles = parts[0]
                if len(parts) == 1:
                    mol_id = f"mol_{lineno:04d}"
                elif len(parts) == 2:
                    mol_id = parts[1]
                else:
                    # 3+ columns: treat the last token as the name
                    # (handles: SMILES index name  or  SMILES name extra…)
                    mol_id = parts[-1]
                molecules.append((smiles, mol_id))
        logger.info("  Parsed %d molecule(s) from file.", len(molecules))

    elif input_p.is_dir():
        logger.info("Batch input: SDF directory %s", input_p)
        sdf_files = sorted(input_p.glob("*.sdf"))
        if not sdf_files:
            logger.warning("No .sdf files found in directory %s", input_p)
        for sdf_file in sdf_files:
            mol_id = sdf_file.stem
            supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=True)
            mol = None
            for candidate in supplier:
                if candidate is not None:
                    mol = candidate
                    break
            if mol is None:
                logger.warning(
                    "Could not read a valid molecule from %s – skipping.",
                    sdf_file.name,
                )
                continue
            smiles = Chem.MolToSmiles(mol)
            molecules.append((smiles, mol_id))
            logger.debug("  %s → %s", mol_id, smiles)
        logger.info("  Collected %d molecule(s) from directory.", len(molecules))

    else:
        raise ValueError(
            f"input_path must be a .txt/.smi file or a directory; got: {input_path!r}"
        )

    # ------------------------------------------------------------------
    # 2. Process each molecule
    # ------------------------------------------------------------------
    os.makedirs(output_dir, exist_ok=True)
    all_dfs: list[pd.DataFrame] = []

    for mol_idx, (smiles, mol_id) in enumerate(molecules, start=1):
        logger.info(
            "━" * 60 + "\nBatch molecule %d/%d: %s (%s)\n" + "━" * 60,
            mol_idx, len(molecules), mol_id, smiles,
        )

        mol_dir = os.path.join(output_dir, mol_id)
        mol_csv = os.path.join(mol_dir, f"{mol_id}_results.csv")

        try:
            df = run_quest_pipeline(
                smiles=smiles,
                receptor_pdb=receptor_pdb,
                receptor_pdbqt=receptor_pdbqt,
                size=size,
                output_csv=mol_csv,
                output_dir=mol_dir,
                **kwargs,
            )
            df.insert(0, "Ligand_ID", mol_id)
            all_dfs.append(df)
            logger.info(
                "Molecule %s finished – %d pose(s) recorded.", mol_id, len(df)
            )

        except Exception as exc:  # noqa: BLE001 – intentionally broad
            logger.error(
                "Pipeline failed for molecule %s (%s): %s",
                mol_id, smiles, exc, exc_info=True,
            )
            continue

    # ------------------------------------------------------------------
    # 3. Assemble master DataFrame
    # ------------------------------------------------------------------
    if not all_dfs:
        logger.warning(
            "Batch run produced no results – every molecule failed or produced "
            "no poses."
        )
        return pd.DataFrame()

    master_df = pd.concat(all_dfs, ignore_index=True)

    output_csv_abs = os.path.abspath(output_csv)
    master_df.to_csv(output_csv_abs, index=False)
    logger.info(
        "Batch complete – %d total pose(s) across %d ligand(s). "
        "Master CSV saved → %s",
        len(master_df), len(all_dfs), output_csv_abs,
    )

    return master_df
