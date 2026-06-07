"""
QUEST – Module 3: The Docking Engine
======================================
Automated pocket detection + ligand docking pipeline using P2Rank and
AutoDock Vina, with OpenBabel handling all PDBQT conversions.

Workflow
--------
1. run_p2rank              – detect binding pockets, return top-pocket centre.
2. prepare_ligand_pdbqt    – SDF → PDBQT via OpenBabel.
3. run_vina                – perform AutoDock Vina docking.
4. convert_vina_output_to_sdf – Vina PDBQT output → multi-pose SDF.
5. dock_ligand             – orchestrator that combines all of the above.

P2Rank is bundled with QUEST at ``tools/p2rank/prank``.  The executable is
resolved in this priority order:

    1. ``QUEST_P2RANK`` environment variable (absolute path to ``prank``).
    2. ``<project_root>/tools/p2rank/prank``  (bundled, works out-of-the-box).
    3. ``prank``  on system PATH (for system-wide installs).
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helper – p2rank executable resolver
# ---------------------------------------------------------------------------

def _resolve_p2rank(hint: str = "prank") -> str:
    """Return the absolute path to the ``prank`` executable.

    Resolution order:
    1. ``QUEST_P2RANK`` environment variable.
    2. ``<quest_package_root>/../tools/p2rank/prank`` (bundled copy).
    3. The value of *hint* (defaults to ``"prank"``), checked on PATH.

    Raises
    ------
    FileNotFoundError
        If none of the above locations yield an executable file.
    """
    # 1 – explicit env override
    env_path = os.environ.get("QUEST_P2RANK", "")
    if env_path and os.path.isfile(env_path):
        logger.debug("P2Rank resolved via QUEST_P2RANK env var: %s", env_path)
        return env_path

    # 2 – bundled copy: quest/docking_engine.py → ../ → tools/p2rank/prank
    bundled = Path(__file__).resolve().parent.parent / "tools" / "p2rank" / "prank"
    if bundled.is_file():
        logger.debug("P2Rank resolved via bundled copy: %s", bundled)
        return str(bundled)

    # 3 – system PATH
    on_path = shutil.which(hint)
    if on_path:
        logger.debug("P2Rank resolved via PATH: %s", on_path)
        return on_path

    raise FileNotFoundError(
        "Cannot locate the P2Rank executable ('prank'). "
        "Either:\n"
        "  a) Set the QUEST_P2RANK environment variable to its absolute path, or\n"
        f" b) Ensure the bundled copy exists at: {bundled}, or\n"
        "  c) Add prank to your system PATH."
    )


# ---------------------------------------------------------------------------
# Function 1
# ---------------------------------------------------------------------------

def run_p2rank(
    pdb_path: str,
    p2rank_exec: str = "prank",
    output_dir: str = ".",
) -> tuple[float, float, float]:
    """Run P2Rank pocket detection and return the top-pocket centre.

    P2Rank is called with ``predict -f {pdb_path} -o {output_dir}/p2rank-pocket``
    and *all* P2Rank outputs (CSV, visualisations, residue scores, etc.) are
    preserved there for inspection.  The CSV is then parsed to extract the
    centre coordinates of the highest-scoring pocket (first data row, ranked
    by score descending).

    Parameters
    ----------
    pdb_path:
        Path to the receptor PDB file.
    p2rank_exec:
        Path or name of the ``prank`` executable.  Resolved automatically
        when the default ``"prank"`` is left unchanged.
    output_dir:
        Parent directory under which a ``p2rank-pocket/`` subdirectory is
        created to hold all P2Rank output files.

    Returns
    -------
    tuple[float, float, float]
        ``(center_x, center_y, center_z)`` of the top-ranked pocket in Å.

    Raises
    ------
    FileNotFoundError
        If the predictions CSV is not generated (P2Rank failed silently).
    RuntimeError
        If P2Rank exits with a non-zero return code.
    ValueError
        If the CSV exists but cannot be parsed.
    """
    exec_path = _resolve_p2rank(p2rank_exec)
    pdb_abs = os.path.abspath(pdb_path)
    pdb_name = os.path.basename(pdb_abs)

    # Create a persistent output directory so all P2Rank files are kept
    p2rank_out_dir = os.path.abspath(os.path.join(output_dir, "p2rank-pocket"))
    os.makedirs(p2rank_out_dir, exist_ok=True)
    logger.info("P2Rank outputs will be saved to: %s", p2rank_out_dir)

    cmd: list[str] = [exec_path, "predict", "-f", pdb_abs, "-o", p2rank_out_dir]

    logger.info("Launching P2Rank: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
            check=True,
        )
        logger.info("P2Rank finished for %s", pdb_name)
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"P2Rank timed out after 600 s for receptor: {pdb_path}"
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"P2Rank exited with code {exc.returncode} for {pdb_path}.\n"
            f"--- stderr ---\n{exc.stderr}\n"
            f"--- stdout ---\n{exc.stdout}"
        ) from exc

    # P2Rank names the CSV after the input PDB: {pdb_name}_predictions.csv
    csv_path = os.path.join(p2rank_out_dir, f"{pdb_name}_predictions.csv")
    if not os.path.isfile(csv_path):
        # Fallback: scan for any *_predictions.csv in case of naming quirks
        matches = [
            f for f in os.listdir(p2rank_out_dir) if f.endswith("_predictions.csv")
        ]
        if not matches:
            raise FileNotFoundError(
                f"P2Rank did not generate a predictions CSV in {p2rank_out_dir}. "
                f"P2Rank stdout:\n" + "\n".join(result.stdout.splitlines()[-30:])
            )
        csv_path = os.path.join(p2rank_out_dir, matches[0])
        logger.debug("Found predictions CSV via fallback scan: %s", csv_path)

    try:
        with open(csv_path, newline="", encoding="utf-8") as fh:
            # P2Rank CSV has leading spaces in header names − strip them
            reader = csv.DictReader(fh, skipinitialspace=True)
            rows = list(reader)

        if not rows:
            raise ValueError(
                f"P2Rank predictions CSV is empty (no pockets detected): {csv_path}"
            )

        top = rows[0]
        center_x = float(top["center_x"])
        center_y = float(top["center_y"])
        center_z = float(top["center_z"])

    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"Could not parse pocket centre from P2Rank CSV {csv_path}: {exc}"
        ) from exc

    logger.info(
        "Top pocket centre: (%.3f, %.3f, %.3f) — full output in %s",
        center_x, center_y, center_z, p2rank_out_dir,
    )
    return (center_x, center_y, center_z)


# ---------------------------------------------------------------------------
# Function 2
# ---------------------------------------------------------------------------

def prepare_ligand_pdbqt(sdf_path: str, pdbqt_path: str) -> None:
    """Convert a ligand SDF to PDBQT format using OpenBabel.

    The SDF is expected to already contain a 3-D conformer (e.g. from
    :mod:`quest.state_generator`).

    Parameters
    ----------
    sdf_path:
        Input SDF file path.
    pdbqt_path:
        Destination PDBQT file path.

    Raises
    ------
    RuntimeError
        If OpenBabel exits with a non-zero return code or produces an empty
        output file.
    """
    cmd: list[str] = [
        "obabel",
        "-isdf", os.path.abspath(sdf_path),
        "-opdbqt",
        "-O", os.path.abspath(pdbqt_path),
        "--partialcharge", "gasteiger",
    ]

    logger.info("Converting ligand SDF → PDBQT: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"OpenBabel SDF→PDBQT conversion failed (exit {exc.returncode}) "
            f"for {sdf_path}.\nstderr: {exc.stderr}"
        ) from exc

    if not os.path.isfile(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
        raise RuntimeError(
            f"OpenBabel produced no output at {pdbqt_path}. "
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

    logger.info("Ligand PDBQT written → %s", pdbqt_path)


# ---------------------------------------------------------------------------
# Function 3
# ---------------------------------------------------------------------------

def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    output_pdbqt: str,
    exhaustiveness: int = 16,
    cwd: str | None = None,
) -> str:
    """Run AutoDock Vina and return the docked poses as a PDBQT string.

    Parameters
    ----------
    receptor_pdbqt:
        Path to the prepared receptor PDBQT.
    ligand_pdbqt:
        Path to the prepared ligand PDBQT.
    center:
        ``(cx, cy, cz)`` – centre of the docking search box in Å.
    size:
        ``(sx, sy, sz)`` – dimensions of the docking search box in Å.
    output_pdbqt:
        Path where Vina writes the docked poses PDBQT.
    exhaustiveness:
        Vina exhaustiveness parameter (higher = more thorough, slower).

    Returns
    -------
    str
        Vina standard output (contains binding affinity table).

    Raises
    ------
    RuntimeError
        If Vina exits with a non-zero return code or times out.
    """
    cx, cy, cz = center
    sx, sy, sz = size

    # Use the absolute path for receptor/ligand (they may live anywhere),
    # but pass ONLY the basename for --out and run Vina with cwd=output_dir.
    # This avoids Vina silently ignoring --out when the absolute path contains
    # spaces (a known Vina argument-parser quirk on Linux).
    out_dir = cwd or os.path.dirname(os.path.abspath(output_pdbqt))
    out_filename = os.path.basename(output_pdbqt)

    cmd: list[str] = [
        "vina",
        "--receptor", os.path.abspath(receptor_pdbqt),
        "--ligand",   os.path.abspath(ligand_pdbqt),
        "--center_x", f"{cx:.4f}",
        "--center_y", f"{cy:.4f}",
        "--center_z", f"{cz:.4f}",
        "--size_x",   f"{sx:.4f}",
        "--size_y",   f"{sy:.4f}",
        "--size_z",   f"{sz:.4f}",
        "--exhaustiveness", str(exhaustiveness),
        "--out", out_filename,
    ]

    logger.info("Launching Vina (cwd=%s): %s", out_dir, " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
            cwd=out_dir,
            check=True,
        )
        final_path = os.path.join(out_dir, out_filename)
        if not os.path.isfile(final_path):
            raise RuntimeError(
                f"Vina exited cleanly but '{out_filename}' was not written to {out_dir}.\n"
                f"stdout:\n{result.stdout}"
            )
        logger.info("Vina completed successfully → %s", final_path)
        return result.stdout

    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"Vina timed out after 3600 s. "
            f"Receptor: {receptor_pdbqt}, Ligand: {ligand_pdbqt}"
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Vina exited with code {exc.returncode}.\n"
            f"--- stderr ---\n{exc.stderr}\n"
            f"--- stdout ---\n{exc.stdout}"
        ) from exc


# ---------------------------------------------------------------------------
# Function 4
# ---------------------------------------------------------------------------

def convert_vina_output_to_sdf(
    vina_out_pdbqt: str,
    output_sdf: str,
) -> None:
    """Convert a Vina multi-pose PDBQT output to SDF using OpenBabel.

    Parameters
    ----------
    vina_out_pdbqt:
        Path to the Vina output PDBQT (may contain multiple poses).
    output_sdf:
        Destination SDF path.

    Raises
    ------
    RuntimeError
        If OpenBabel fails or produces an empty SDF.
    """
    cmd: list[str] = [
        "obabel",
        "-ipdbqt", os.path.abspath(vina_out_pdbqt),
        "-osdf",
        "-O", os.path.abspath(output_sdf),
    ]

    logger.info("Converting Vina PDBQT → SDF: %s", " ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"OpenBabel PDBQT→SDF conversion failed (exit {exc.returncode}) "
            f"for {vina_out_pdbqt}.\nstderr: {exc.stderr}"
        ) from exc

    if not os.path.isfile(output_sdf) or os.path.getsize(output_sdf) == 0:
        raise RuntimeError(
            f"OpenBabel produced no SDF output at {output_sdf}. "
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

    logger.info("Docked SDF written → %s", output_sdf)


# ---------------------------------------------------------------------------
# Function 5
# ---------------------------------------------------------------------------

def dock_ligand(
    qm_sdf: str,
    receptor_pdbqt: str,
    size: tuple[float, float, float],
    center: tuple[float, float, float] | None = None,
    p2rank_pdb: str | None = None,
    p2rank_exec: str = "prank",
    output_dir: str = ".",
    exhaustiveness: int = 16,
) -> str:
    """Orchestrate the full docking pipeline for a single ligand SDF.

    Either *center* must be provided directly, or *p2rank_pdb* must be given
    so that P2Rank can detect the binding pocket automatically.

    Pipeline
    --------
    .. code-block:: text

        qm_sdf
          │
          ├─ [if center is None] run_p2rank(p2rank_pdb) → center
          │
          ├─ prepare_ligand_pdbqt()  → ligand.pdbqt (temp)
          ├─ run_vina()              → docked.pdbqt
          ├─ convert_vina_output_to_sdf() → docked_poses.sdf
          └─ return absolute path to docked_poses.sdf

    Parameters
    ----------
    qm_sdf:
        Path to the QM-optimised ligand SDF (output of Module 2).
    receptor_pdbqt:
        Path to the prepared receptor PDBQT file.
    size:
        ``(sx, sy, sz)`` – docking box dimensions in Å.
        A value of ``(20.0, 20.0, 20.0)`` is a reasonable starting point.
    center:
        ``(cx, cy, cz)`` – docking box centre in Å.  Pass ``None`` to have
        P2Rank determine the centre automatically.
    p2rank_pdb:
        Path to the receptor PDB used by P2Rank.  Required when *center*
        is ``None``.
    p2rank_exec:
        Path or name of the ``prank`` executable (auto-resolved by default).
    output_dir:
        Directory where ``docked.pdbqt`` and ``docked_poses.sdf`` are saved.
    exhaustiveness:
        Vina exhaustiveness (default 16 – good balance of speed/accuracy).

    Returns
    -------
    str
        Absolute path to ``docked_poses.sdf``.

    Raises
    ------
    ValueError
        If both *center* and *p2rank_pdb* are ``None``.
    RuntimeError
        If any step in the pipeline fails.
    """
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1 – resolve box centre
    # ------------------------------------------------------------------
    if center is None:
        if p2rank_pdb is None:
            raise ValueError(
                "Either 'center' coordinates or 'p2rank_pdb' must be provided. "
                "Cannot determine the docking box centre without one of these."
            )
        logger.info("Running P2Rank to detect binding pocket …")
        center = run_p2rank(p2rank_pdb, p2rank_exec, output_dir=output_dir)
        logger.info("P2Rank pocket centre: %s", center)

    # ------------------------------------------------------------------
    # Step 2 – ligand SDF → PDBQT (temporary file)
    # ------------------------------------------------------------------
    stem = Path(qm_sdf).stem
    tmp_ligand_pdbqt = os.path.abspath(
        os.path.join(output_dir, f"_{stem}_ligand_tmp.pdbqt")
    )
    prepare_ligand_pdbqt(qm_sdf, tmp_ligand_pdbqt)

    # ------------------------------------------------------------------
    # Step 3 – Vina docking
    # ------------------------------------------------------------------
    docked_pdbqt = os.path.abspath(os.path.join(output_dir, f"{stem}_docked.pdbqt"))

    try:
        vina_stdout = run_vina(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=tmp_ligand_pdbqt,
            center=center,
            size=size,
            output_pdbqt=docked_pdbqt,
            exhaustiveness=exhaustiveness,
            cwd=output_dir,
        )
        # Log Vina's affinity table
        for line in vina_stdout.splitlines():
            logger.info("[Vina] %s", line)
    finally:
        # Always clean up the temporary ligand PDBQT
        try:
            os.unlink(tmp_ligand_pdbqt)
        except OSError:
            pass

    # ------------------------------------------------------------------
    # Step 4 – Vina PDBQT → SDF
    # ------------------------------------------------------------------
    docked_sdf = os.path.abspath(os.path.join(output_dir, f"{stem}_docked_poses.sdf"))
    convert_vina_output_to_sdf(docked_pdbqt, docked_sdf)

    logger.info("dock_ligand complete → %s", docked_sdf)
    return docked_sdf
