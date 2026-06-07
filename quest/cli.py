"""
QUEST command-line interface
============================
Entry point registered as the ``quest`` console script.

Usage examples
--------------
  quest strain ref_ligand.sdf docked_pose.sdf
  quest states "CC(=O)Nc1ccc(O)cc1" -n 3 -T 8 -o ./states
  quest dock tautomer_0_qm.sdf receptor.pdbqt --p2rank-pdb receptor.pdb
  quest run "CC(=O)Nc1ccc(O)cc1" receptor.pdb receptor.pdbqt --size 20 20 20
"""

from __future__ import annotations

import logging
import os
import sys

import click

from quest import __version__


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(version=__version__, prog_name="QUEST")
@click.option(
    "--log-level",
    default="WARNING",
    show_default=True,
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False
    ),
    help="Set the logging verbosity.",
)
def main(log_level: str) -> None:
    """QUEST – Quantum Energetics Strain Tool.

    A computational chemistry CLI for evaluating conformational strain
    energies of docked ligand poses using GFN2-xTB.
    """
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s  %(name)-30s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# `quest strain` sub-command
# ---------------------------------------------------------------------------

@main.command("strain")
@click.argument("ref_sdf", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument(
    "docked_sdf", type=click.Path(exists=True, dir_okay=False, readable=True)
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Print intermediate energies and xTB command details.",
)
@click.option(
    "--output-dir",
    "-o",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, writable=True),
    help="Directory where temporary files (docked_ready.xyz, constrain.inp, "
         "ref_free.xyz) are written.",
)
def strain_cmd(
    ref_sdf: str,
    docked_sdf: str,
    verbose: bool,
    output_dir: str,
) -> None:
    """Calculate the conformational strain energy of a docked ligand pose.

    \b
    Arguments
    ---------
    REF_SDF    Path to the reference (free) ligand SDF with a 3-D conformer.
    DOCKED_SDF Path to the docked pose SDF from the docking engine.

    \b
    Output
    ------
    Prints the strain energy in kcal mol⁻¹ to stdout.
    A positive value means the docked pose is strained relative to the
    relaxed free ligand.

    \b
    Examples
    --------
      quest strain ref.sdf docked.sdf
      quest strain ref.sdf docked.sdf -v --output-dir ./scratch
    """
    # --- promote to INFO when --verbose is set --------------------------------
    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    # --- lazy import so the CLI launches fast even without rdkit installed ----
    try:
        from quest.strain_evaluator import calculate_strain
    except ImportError as exc:
        click.echo(
            f"[QUEST ERROR] Could not import strain_evaluator: {exc}\n"
            f"Make sure rdkit and xtb are available in your environment.",
            err=True,
        )
        sys.exit(1)

    # --- ensure output directory exists --------------------------------------
    os.makedirs(output_dir, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(output_dir)

    # Convert paths to absolute *before* chdir so relative inputs still work.
    ref_sdf_abs = os.path.abspath(os.path.join(original_dir, ref_sdf))
    docked_sdf_abs = os.path.abspath(os.path.join(original_dir, docked_sdf))

    ref_name = os.path.basename(ref_sdf_abs)
    docked_name = os.path.basename(docked_sdf_abs)

    click.echo(f"QUEST  v{__version__}  –  Strain Evaluator")
    click.echo(f"  Reference : {ref_name}")
    click.echo(f"  Docked    : {docked_name}")
    click.echo(f"  Output dir: {os.path.abspath(output_dir)}")
    click.echo("")

    try:
        strain_kcal: float = calculate_strain(ref_sdf_abs, docked_sdf_abs)
    except (ValueError, IOError, RuntimeError) as exc:
        click.echo(f"[QUEST ERROR] {exc}", err=True)
        os.chdir(original_dir)
        sys.exit(1)
    finally:
        os.chdir(original_dir)

    click.echo(f"  Strain energy: {strain_kcal:+.4f} kcal mol⁻¹")

    if verbose:
        click.echo("")
        if strain_kcal < 1.5:
            click.echo("  Interpretation: low strain (<1.5 kcal mol⁻¹)  ✓")
        elif strain_kcal < 5.0:
            click.echo("  Interpretation: moderate strain (1.5–5.0 kcal mol⁻¹)  ⚠")
        else:
            click.echo("  Interpretation: high strain (>5.0 kcal mol⁻¹)  ✗")


# ---------------------------------------------------------------------------
# `quest states` sub-command
# ---------------------------------------------------------------------------

@main.command("states")
@click.argument("smiles")
@click.option(
    "-n",
    "--max-tautomers",
    default=3,
    show_default=True,
    type=click.IntRange(1, 20),
    help="Maximum number of tautomers to enumerate and process.",
)
@click.option(
    "-T",
    "--threads",
    default=4,
    show_default=True,
    type=click.IntRange(1, 256),
    help="CPU threads passed to CREST via -T.",
)
@click.option(
    "-o",
    "--output-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, writable=True),
    help="Directory where tautomer_{i}_qm.sdf files are written.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Promote logging to INFO so CREST progress is visible.",
)
def states_cmd(
    smiles: str,
    max_tautomers: int,
    threads: int,
    output_dir: str,
    verbose: bool,
) -> None:
    """Enumerate tautomers and find QM minimum geometries via CREST.

    \b
    Arguments
    ---------
    SMILES   Input SMILES string of the ligand (quote if it contains brackets).

    \b
    Output
    ------
    Writes one SDF per tautomer: tautomer_0_qm.sdf, tautomer_1_qm.sdf, …
    Prints the list of saved paths to stdout.

    \b
    Examples
    --------
      quest states "CC(=O)Nc1ccc(O)cc1"
      quest states "CC(=O)Nc1ccc(O)cc1" -n 5 -T 8 -o ./states_out -v
    """
    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    try:
        from quest.state_generator import prep_ligand_states
    except ImportError as exc:
        click.echo(
            f"[QUEST ERROR] Could not import state_generator: {exc}\n"
            f"Make sure rdkit and crest are available in your environment.",
            err=True,
        )
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    output_dir_abs = os.path.abspath(output_dir)

    click.echo(f"QUEST  v{__version__}  –  State Generator")
    click.echo(f"  SMILES      : {smiles}")
    click.echo(f"  Max tautomers: {max_tautomers}")
    click.echo(f"  Threads     : {threads}")
    click.echo(f"  Output dir  : {output_dir_abs}")
    click.echo("")

    try:
        saved: list[str] = prep_ligand_states(
            smiles,
            output_dir=output_dir_abs,
            max_tautomers=max_tautomers,
            threads=threads,
        )
    except (ValueError, RuntimeError) as exc:
        click.echo(f"[QUEST ERROR] {exc}", err=True)
        sys.exit(1)

    click.echo(f"  Completed – {len(saved)} SDF file(s) written:")
    for path in saved:
        click.echo(f"    {path}")


# ---------------------------------------------------------------------------
# `quest dock` sub-command
# ---------------------------------------------------------------------------

@main.command("dock")
@click.argument(
    "qm_sdf",
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.argument(
    "receptor_pdbqt",
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.option(
    "--size",
    nargs=3,
    type=float,
    default=(20.0, 20.0, 20.0),
    show_default=True,
    metavar="SX SY SZ",
    help="Docking box dimensions in Å (x y z).",
)
@click.option(
    "--center",
    nargs=3,
    type=float,
    default=None,
    metavar="CX CY CZ",
    help="Docking box centre in Å (x y z). If omitted, P2Rank auto-detects it.",
)
@click.option(
    "--p2rank-pdb",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Receptor PDB for automatic pocket detection via P2Rank (required when --center is omitted).",
)
@click.option(
    "--exhaustiveness",
    "-e",
    default=16,
    show_default=True,
    type=click.IntRange(1, 128),
    help="Vina exhaustiveness (higher = more thorough).",
)
@click.option(
    "-o",
    "--output-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, writable=True),
    help="Directory where docked PDBQT and SDF are written.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Promote logging to INFO so Vina/P2Rank progress is visible.",
)
def dock_cmd(
    qm_sdf: str,
    receptor_pdbqt: str,
    size: tuple,
    center: tuple | None,
    p2rank_pdb: str | None,
    exhaustiveness: int,
    output_dir: str,
    verbose: bool,
) -> None:
    """Dock a ligand SDF into a prepared receptor using AutoDock Vina.

    \b
    Arguments
    ---------
    QM_SDF          QM-optimised ligand SDF (output of 'quest states').
    RECEPTOR_PDBQT  Prepared receptor PDBQT file.

    \b
    Centre resolution (in order of priority)
    -----------------------------------------
      1. --center CX CY CZ   (explicit coordinates)
      2. --p2rank-pdb PDB    (automatic pocket detection via P2Rank)

    \b
    Output
    ------
    Writes {stem}_docked.pdbqt and {stem}_docked_poses.sdf to --output-dir.
    Prints the final SDF path to stdout.

    \b
    Examples
    --------
      quest dock tautomer_0_qm.sdf receptor.pdbqt --p2rank-pdb receptor.pdb
      quest dock tautomer_0_qm.sdf receptor.pdbqt --center 10.5 -3.2 22.0 --size 22 22 22
      quest dock tautomer_0_qm.sdf receptor.pdbqt --p2rank-pdb receptor.pdb -e 32 -v -o ./dock_out
    """
    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    try:
        from quest.docking_engine import dock_ligand
    except ImportError as exc:
        click.echo(
            f"[QUEST ERROR] Could not import docking_engine: {exc}\n"
            "Make sure openbabel and vina are available in your environment.",
            err=True,
        )
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    output_dir_abs = os.path.abspath(output_dir)

    center_tuple: tuple[float, float, float] | None = (
        (float(center[0]), float(center[1]), float(center[2])) if center else None
    )
    size_tuple: tuple[float, float, float] = (
        float(size[0]), float(size[1]), float(size[2])
    )

    click.echo(f"QUEST  v{__version__}  –  Docking Engine")
    click.echo(f"  Ligand SDF      : {qm_sdf}")
    click.echo(f"  Receptor PDBQT  : {receptor_pdbqt}")
    click.echo(f"  Box size        : {size_tuple}")
    click.echo(
        f"  Box centre      : "
        + (f"{center_tuple}" if center_tuple else "auto (P2Rank)")
    )
    click.echo(f"  Exhaustiveness  : {exhaustiveness}")
    click.echo(f"  Output dir      : {output_dir_abs}")
    click.echo("")

    try:
        docked_sdf = dock_ligand(
            qm_sdf=os.path.abspath(qm_sdf),
            receptor_pdbqt=os.path.abspath(receptor_pdbqt),
            size=size_tuple,
            center=center_tuple,
            p2rank_pdb=os.path.abspath(p2rank_pdb) if p2rank_pdb else None,
            output_dir=output_dir_abs,
            exhaustiveness=exhaustiveness,
        )
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        click.echo(f"[QUEST ERROR] {exc}", err=True)
        sys.exit(1)

    click.echo(f"  Docking complete → {docked_sdf}")


# ---------------------------------------------------------------------------
# `quest run` sub-command  (master orchestrator)
# ---------------------------------------------------------------------------

@main.command("run")
@click.argument("smiles")
@click.argument(
    "receptor_pdbqt",
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.option(
    "--p2rank-pdb",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Receptor PDB fed to P2Rank for automatic binding-pocket detection.",
)
@click.option(
    "--size",
    nargs=3,
    type=float,
    default=(20.0, 20.0, 20.0),
    show_default=True,
    metavar="SX SY SZ",
    help="Vina docking box dimensions in Å.",
)
@click.option(
    "-n",
    "--max-tautomers",
    default=3,
    show_default=True,
    type=click.IntRange(1, 20),
    help="Maximum tautomers to enumerate.",
)
@click.option(
    "-T",
    "--threads",
    default=4,
    show_default=True,
    type=click.IntRange(1, 256),
    help="CPU threads for CREST.",
)
@click.option(
    "-e",
    "--exhaustiveness",
    default=16,
    show_default=True,
    type=click.IntRange(1, 128),
    help="Vina exhaustiveness.",
)
@click.option(
    "--strain-cutoff",
    default=7.0,
    show_default=True,
    type=float,
    help="Strain threshold (kcal mol⁻¹) for the Viable flag.",
)
@click.option(
    "--output-csv",
    default="quest_results.csv",
    show_default=True,
    help="Path for the output results CSV.",
)
@click.option(
    "-o",
    "--output-dir",
    default="quest_output",
    show_default=True,
    type=click.Path(file_okay=False, writable=True),
    help="Root directory for all pipeline output files.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Promote logging to INFO for full pipeline progress.",
)
def run_cmd(
    smiles: str,
    receptor_pdbqt: str,
    p2rank_pdb: str,
    size: tuple,
    max_tautomers: int,
    threads: int,
    exhaustiveness: int,
    strain_cutoff: float,
    output_csv: str,
    output_dir: str,
    verbose: bool,
) -> None:
    """Run the full QUEST pipeline: states → dock → strain → CSV.

    \b
    Arguments
    ---------
    SMILES          Ligand SMILES string (quote if it contains brackets).
    RECEPTOR_PDBQT  Prepared receptor PDBQT passed to AutoDock Vina.

    \b
    Pipeline stages
    ---------------
    1. CREST/GFN2  – enumerate tautomers, QM-minimise each (Module 2)
    2. P2Rank      – detect top binding pocket from --p2rank-pdb
    3. Vina        – dock each tautomer into the pocket (Module 3)
    4. xTB GFN2    – calculate conformational strain per pose (Module 1)

    \b
    Output
    ------
    quest_results.csv (or --output-csv path) with columns:
      Tautomer | Pose | Vina_Affinity | QM_Strain | Viable
    Sorted by Vina_Affinity ascending (best first).

    \b
    Examples
    --------
      quest run "CC(=O)Oc1ccccc1C(=O)O" receptor.pdbqt --p2rank-pdb receptor.pdb --size 20 20 20
      quest run "CC(=O)Oc1ccccc1C(=O)O" receptor.pdbqt --p2rank-pdb receptor.pdb \\
          --size 22 22 22 -n 5 -T 12 -e 32 --strain-cutoff 5.0 -v -o ./run1
    """
    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    try:
        from quest.pipeline import run_quest_pipeline
    except ImportError as exc:
        click.echo(
            f"[QUEST ERROR] Could not import pipeline: {exc}\n"
            "Make sure rdkit, crest, xtb, vina, and openbabel are installed.",
            err=True,
        )
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    output_dir_abs = os.path.abspath(output_dir)
    output_csv_abs = (
        os.path.abspath(output_csv)
        if os.path.isabs(output_csv)
        else os.path.abspath(os.path.join(output_dir_abs, output_csv))
    )
    size_tuple = (float(size[0]), float(size[1]), float(size[2]))

    click.echo(f"QUEST  v{__version__}  –  Full Pipeline")
    click.echo(f"  SMILES          : {smiles}")
    click.echo(f"  P2Rank PDB      : {p2rank_pdb}  (pocket detection)")
    click.echo(f"  Receptor PDBQT  : {receptor_pdbqt}  (Vina docking)")
    click.echo(f"  Box size        : {size_tuple}")
    click.echo(f"  Max tautomers   : {max_tautomers}")
    click.echo(f"  CREST threads   : {threads}")
    click.echo(f"  Exhaustiveness  : {exhaustiveness}")
    click.echo(f"  Strain cutoff   : {strain_cutoff} kcal mol⁻¹")
    click.echo(f"  Output dir      : {output_dir_abs}")
    click.echo(f"  Results CSV     : {output_csv_abs}")
    click.echo("")

    try:
        df = run_quest_pipeline(
            smiles=smiles,
            receptor_pdb=os.path.abspath(p2rank_pdb),
            receptor_pdbqt=os.path.abspath(receptor_pdbqt),
            size=size_tuple,
            output_csv=output_csv_abs,
            output_dir=output_dir_abs,
            strain_cutoff=strain_cutoff,
            max_tautomers=max_tautomers,
            threads=threads,
            exhaustiveness=exhaustiveness,
        )
    except (ValueError, RuntimeError) as exc:
        click.echo(f"[QUEST ERROR] {exc}", err=True)
        sys.exit(1)

    click.echo(f"\n  Pipeline complete – {len(df)} pose(s) scored.")
    click.echo(f"  Results saved → {output_csv_abs}")
    click.echo("")

    # Pretty-print summary table
    viable = df[df["Viable"] == True]  # noqa: E712
    click.echo(f"  Viable poses (strain ≤ {strain_cutoff} kcal mol⁻¹): {len(viable)}")
    if not df.empty:
        click.echo("")
        click.echo(df.to_string(index=True))


# ---------------------------------------------------------------------------
# `quest batch` sub-command
# ---------------------------------------------------------------------------

@main.command("batch")
@click.argument(
    "input_path",
    type=click.Path(exists=True, readable=True),
)
@click.argument(
    "receptor_pdbqt",
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.option(
    "--p2rank-pdb",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Receptor PDB fed to P2Rank for automatic binding-pocket detection.",
)
@click.option(
    "--size",
    nargs=3,
    type=float,
    default=(20.0, 20.0, 20.0),
    show_default=True,
    metavar="SX SY SZ",
    help="Vina docking box dimensions in Å.",
)
@click.option(
    "-n",
    "--max-tautomers",
    default=3,
    show_default=True,
    type=click.IntRange(1, 20),
    help="Maximum tautomers to enumerate per ligand.",
)
@click.option(
    "-T",
    "--threads",
    default=4,
    show_default=True,
    type=click.IntRange(1, 256),
    help="CPU threads for CREST.",
)
@click.option(
    "-e",
    "--exhaustiveness",
    default=16,
    show_default=True,
    type=click.IntRange(1, 128),
    help="Vina exhaustiveness.",
)
@click.option(
    "--strain-cutoff",
    default=7.0,
    show_default=True,
    type=float,
    help="Strain threshold (kcal mol⁻¹) for the Viable flag.",
)
@click.option(
    "--output-csv",
    default="quest_batch_results.csv",
    show_default=True,
    help="Path for the master batch results CSV.",
)
@click.option(
    "-o",
    "--output-dir",
    default="quest_batch_output",
    show_default=True,
    type=click.Path(file_okay=False, writable=True),
    help="Root directory for all batch output files. A sub-directory is "
         "created per ligand.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Promote logging to INFO for full pipeline progress.",
)
def batch_cmd(
    input_path: str,
    receptor_pdbqt: str,
    p2rank_pdb: str,
    size: tuple,
    max_tautomers: int,
    threads: int,
    exhaustiveness: int,
    strain_cutoff: float,
    output_csv: str,
    output_dir: str,
    verbose: bool,
) -> None:
    """Run the full QUEST pipeline over a batch of ligands.

    \b
    Arguments
    ---------
    INPUT_PATH      A .smi/.txt file (one SMILES [ID] per line) OR a directory
                    of .sdf files.
    RECEPTOR_PDBQT  Prepared receptor PDBQT passed to AutoDock Vina.

    \b
    Output
    ------
    One sub-directory per ligand under --output-dir, plus a master
    quest_batch_results.csv (or --output-csv path) combining all poses with
    an extra Ligand_ID column.

    \b
    Examples
    --------
      quest batch ligands.smi receptor.pdbqt --p2rank-pdb receptor.pdb
      quest batch ligands.smi receptor.pdbqt --p2rank-pdb receptor.pdb \\
          --size 22 22 22 -n 2 -T 12 -e 32 -v -o ./batch_run
      quest batch ./sdf_library/ receptor.pdbqt --p2rank-pdb receptor.pdb -v
    """
    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    try:
        from quest.pipeline import run_quest_batch
    except ImportError as exc:
        click.echo(
            f"[QUEST ERROR] Could not import pipeline: {exc}\n"
            "Make sure rdkit, crest, xtb, vina, and openbabel are installed.",
            err=True,
        )
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    output_dir_abs = os.path.abspath(output_dir)
    output_csv_abs = (
        os.path.abspath(output_csv)
        if os.path.isabs(output_csv)
        else os.path.abspath(os.path.join(output_dir_abs, output_csv))
    )
    size_tuple = (float(size[0]), float(size[1]), float(size[2]))

    click.echo(f"QUEST  v{__version__}  –  Batch Pipeline")
    click.echo(f"  Input           : {input_path}")
    click.echo(f"  P2Rank PDB      : {p2rank_pdb}  (pocket detection)")
    click.echo(f"  Receptor PDBQT  : {receptor_pdbqt}  (Vina docking)")
    click.echo(f"  Box size        : {size_tuple}")
    click.echo(f"  Max tautomers   : {max_tautomers}")
    click.echo(f"  CREST threads   : {threads}")
    click.echo(f"  Exhaustiveness  : {exhaustiveness}")
    click.echo(f"  Strain cutoff   : {strain_cutoff} kcal mol⁻¹")
    click.echo(f"  Output dir      : {output_dir_abs}")
    click.echo(f"  Master CSV      : {output_csv_abs}")
    click.echo("")

    try:
        df = run_quest_batch(
            input_path=os.path.abspath(input_path),
            receptor_pdb=os.path.abspath(p2rank_pdb),
            receptor_pdbqt=os.path.abspath(receptor_pdbqt),
            size=size_tuple,
            output_csv=output_csv_abs,
            output_dir=output_dir_abs,
            strain_cutoff=strain_cutoff,
            max_tautomers=max_tautomers,
            threads=threads,
            exhaustiveness=exhaustiveness,
        )
    except (ValueError, RuntimeError) as exc:
        click.echo(f"[QUEST ERROR] {exc}", err=True)
        sys.exit(1)

    if df.empty:
        click.echo("  No poses were scored – check logs for errors.")
        sys.exit(1)

    n_ligands = df["Ligand_ID"].nunique() if "Ligand_ID" in df.columns else "?"
    viable = df[df["Viable"] == True]  # noqa: E712
    click.echo(f"\n  Batch complete – {n_ligands} ligand(s), {len(df)} total pose(s).")
    click.echo(f"  Viable poses (strain ≤ {strain_cutoff} kcal mol⁻¹): {len(viable)}")
    click.echo(f"  Master CSV saved → {output_csv_abs}")
    click.echo("")
    click.echo(df.to_string(index=True))
