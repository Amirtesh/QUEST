"""
QUEST – Module 1: The Strain Evaluator
========================================
Calculates the conformational strain energy of a docked ligand pose
relative to its lowest-energy free conformation using xTB (GFN2-xTB).

Workflow
--------
1. transfer_topology   – transplant docked 3-D coordinates onto the
                         reference heavy-atom graph and re-add H atoms.
2. create_xtb_constraint – write a constrained-optimisation input file
                           that freezes all heavy atoms.
3. run_xtb_gfn2        – subprocess wrapper around the xTB binary with
                         full error-handling and a 20-minute safety
                         timeout.
4. extract_energy      – parse TOTAL ENERGY from xTB stdout.
5. calculate_strain    – orchestrator that combines the above steps and
                         returns the strain energy in kcal mol⁻¹.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile

from rdkit import Chem
from rdkit.Chem import rdFMCS

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Function 1
# ---------------------------------------------------------------------------

def transfer_topology(
    ref_sdf_path: str,
    docked_sdf_path: str,
    output_xyz_path: str,
) -> int:
    """Transfer docked 3-D coordinates onto the reference molecular graph.

    The reference SDF supplies the correct heavy-atom connectivity and
    protonation state; the docked SDF supplies the binding-pose geometry.
    An MCS match is used to establish the atom mapping, after which the
    docked Cartesian coordinates are written back onto the reference
    conformer.  Hydrogens are then re-added with 3-D coordinates and the
    result is written as an XYZ file.

    Parameters
    ----------
    ref_sdf_path:
        Path to the reference (free) ligand SDF – must contain H atoms.
    docked_sdf_path:
        Path to the docked ligand SDF as it comes from the docking engine.
    output_xyz_path:
        Destination for the topology-corrected XYZ file.

    Returns
    -------
    int
        Total number of atoms (including H) in the output XYZ file.

    Raises
    ------
    ValueError
        If the MCS graph-match does not cover all heavy atoms of the
        reference molecule (incomplete match → coordinates would be
        meaningless).
    IOError / RuntimeError
        On any file-read or RDKit failure.
    """
    try:
        # --- load molecules ---------------------------------------------------
        ref_with_H: Chem.Mol = next(
            Chem.SDMolSupplier(ref_sdf_path, removeHs=False)
        )
        if ref_with_H is None:
            raise IOError(f"Could not read reference SDF: {ref_sdf_path}")

        ref: Chem.Mol = Chem.RemoveHs(ref_with_H)

        docked: Chem.Mol = next(
            Chem.SDMolSupplier(docked_sdf_path, sanitize=False)
        )
        if docked is None:
            raise IOError(f"Could not read docked SDF: {docked_sdf_path}")

        # --- maximum common substructure match --------------------------------
        mcs = rdFMCS.FindMCS(
            [ref, docked],
            bondCompare=rdFMCS.BondCompare.CompareAny,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
        )
        patt: Chem.Mol = Chem.MolFromSmarts(mcs.smartsString)

        ref_match: tuple[int, ...] = ref.GetSubstructMatch(patt)
        docked_match: tuple[int, ...] = docked.GetSubstructMatch(patt)

        if len(ref_match) != ref.GetNumAtoms():
            raise ValueError(
                f"Graph match failed. Reference has {ref.GetNumAtoms()} heavy "
                f"atoms but MCS matched only {len(ref_match)}. "
                f"Ensure both SDF files represent the same scaffold."
            )

        # --- overwrite reference conformer with docked coordinates -----------
        ref_conf: Chem.Conformer = ref.GetConformer()
        docked_conf: Chem.Conformer = docked.GetConformer()

        for r_idx, d_idx in zip(ref_match, docked_match):
            ref_conf.SetAtomPosition(r_idx, docked_conf.GetAtomPosition(d_idx))

        # --- re-add H with 3-D coordinates and write XYZ ---------------------
        final_mol: Chem.Mol = Chem.AddHs(ref, addCoords=True)
        Chem.MolToXYZFile(final_mol, output_xyz_path)

        n_atoms: int = final_mol.GetNumAtoms()
        logger.info(
            "Topology transfer complete → %s  (%d atoms)", output_xyz_path, n_atoms
        )
        return n_atoms

    except (IOError, ValueError):
        raise
    except Exception as exc:
        raise RuntimeError(
            f"Unexpected error during topology transfer: {exc}"
        ) from exc


# ---------------------------------------------------------------------------
# Function 2
# ---------------------------------------------------------------------------

def create_xtb_constraint(filepath: str = "constrain.inp") -> None:
    """Write an xTB constraint file that freezes all heavy atoms during optimisation.

    The constraint applies a harmonic force constant of 1.0 a.u. to all
    atoms of the listed elements, effectively holding the heavy-atom
    skeleton in the docked pose while allowing H atoms to relax.

    Parameters
    ----------
    filepath:
        Destination path for the constraint file (default: ``constrain.inp``).

    Raises
    ------
    IOError
        If the file cannot be written.
    """
    content: str = (
        "$constrain\n"
        "  force constant=1.0\n"
        "  elements: C,O,N,S,P,F,Cl,Br,I\n"
        "$end\n"
    )
    try:
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(content)
        logger.info("xTB constraint file written → %s", filepath)
    except OSError as exc:
        raise IOError(f"Could not write constraint file to {filepath}: {exc}") from exc


# ---------------------------------------------------------------------------
# Function 3
# ---------------------------------------------------------------------------

def run_xtb_gfn2(
    xyz_file: str,
    constraint_file: str | None = None,
) -> str:
    """Execute an xTB GFN2-xTB geometry optimisation in implicit water.

    The calculation is run inside a temporary scratch directory so that xTB
    output files do not pollute the working directory.  Both ``stdout`` and
    ``stderr`` are captured.

    Parameters
    ----------
    xyz_file:
        Path to the input XYZ (or SDF) file – converted to an absolute path
        before being passed to xTB.
    constraint_file:
        Optional path to an xTB ``--input`` file (e.g. ``constrain.inp``).
        When provided the constrained optimisation flag is appended to the
        command.

    Returns
    -------
    str
        The full standard-output string produced by xTB (contains energies,
        gradient norms, etc.).

    Raises
    ------
    RuntimeError
        If xTB exits with a non-zero return code **or** does not complete
        within 1200 seconds (20 minutes).
    """
    xyz_abs: str = os.path.abspath(xyz_file)

    cmd: list[str] = [
        "xtb", xyz_abs,
        "--gfn", "2",
        "--alpb", "water",
        "--opt",
    ]

    if constraint_file is not None:
        constraint_abs: str = os.path.abspath(constraint_file)
        cmd.extend(["--input", constraint_abs])

    logger.info("Launching xTB: %s", " ".join(cmd))

    with tempfile.TemporaryDirectory() as scratch_dir:
        try:
            result: subprocess.CompletedProcess = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1200,
                cwd=scratch_dir,
                check=True,
            )
            logger.info(
                "xTB finished successfully for %s", os.path.basename(xyz_file)
            )
            return result.stdout

        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"xTB calculation timed out after 1200 s for input file: {xyz_file}. "
                f"The system may be too large or the geometry is highly unfavourable."
            )

        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"xTB returned non-zero exit code {exc.returncode} "
                f"for input file: {xyz_file}.\n"
                f"--- xTB stderr ---\n{exc.stderr}\n"
                f"--- xTB stdout (last 40 lines) ---\n"
                + "\n".join(exc.stdout.splitlines()[-40:])
            ) from exc


# ---------------------------------------------------------------------------
# Function 4
# ---------------------------------------------------------------------------

def extract_energy(xtb_stdout: str) -> float:
    """Parse the GFN2-xTB total electronic energy from captured stdout.

    xTB prints a summary table at the end of each job.  The relevant line
    has the form::

        | TOTAL ENERGY             -42.123456789012 Eh   |

    This function locates that line, extracts the floating-point value, and
    returns it in Hartrees.

    Parameters
    ----------
    xtb_stdout:
        The complete standard-output string returned by :func:`run_xtb_gfn2`.

    Returns
    -------
    float
        Total electronic energy in Hartrees (Eh).

    Raises
    ------
    ValueError
        If no ``| TOTAL ENERGY`` line is found (e.g. the job crashed before
        printing the summary).
    """
    for line in xtb_stdout.splitlines():
        if "| TOTAL ENERGY" in line:
            # Typical format: | TOTAL ENERGY   -XX.XXXXXXXXXXXX Eh   |
            parts = line.split()
            for part in parts:
                try:
                    return float(part)
                except ValueError:
                    continue

    raise ValueError(
        "Could not locate '| TOTAL ENERGY' in xTB output. "
        "The calculation may not have converged or the output was truncated."
    )


# ---------------------------------------------------------------------------
# Function 5
# ---------------------------------------------------------------------------

def calculate_strain(ref_sdf: str, docked_sdf: str) -> float:
    """Calculate the conformational strain energy of a docked ligand pose.

    The strain energy is defined as:

    .. math::

        \\Delta E_{\\text{strain}} = (E_{\\text{docked}} - E_{\\text{free}})
                                     \\times 627.509 \\; \\text{kcal mol}^{-1}

    where

    * :math:`E_{\\text{docked}}` is the GFN2-xTB energy of the ligand
      locked in its docked geometry (heavy atoms frozen, H atoms
      relaxed),
    * :math:`E_{\\text{free}}` is the GFN2-xTB energy of the same
      ligand after unconstrained gas-phase optimisation from the
      reference geometry.

    Parameters
    ----------
    ref_sdf:
        Path to the reference (free) ligand SDF.  Must contain a 3-D
        conformer that serves as the starting point for the free-energy
        optimisation.
    docked_sdf:
        Path to the docked ligand SDF as produced by the docking engine.

    Returns
    -------
    float
        Strain energy in kcal mol⁻¹.  A positive value indicates that the
        docked pose is higher in energy than the relaxed reference.

    Raises
    ------
    RuntimeError
        If any xTB calculation fails or times out.
    ValueError
        If energies cannot be parsed from xTB output.
    IOError
        If RDKit cannot read the input SDF files.
    """
    # ------------------------------------------------------------------
    # Step 1 – transfer docked coordinates onto the reference topology
    # ------------------------------------------------------------------
    docked_xyz: str = "docked_ready.xyz"
    n_atoms: int = transfer_topology(ref_sdf, docked_sdf, docked_xyz)
    logger.info("Topology transfer produced %d atoms → %s", n_atoms, docked_xyz)

    # ------------------------------------------------------------------
    # Step 2 – write the heavy-atom constraint file
    # ------------------------------------------------------------------
    constraint_file: str = "constrain.inp"
    create_xtb_constraint(constraint_file)

    # ------------------------------------------------------------------
    # Step 3 – constrained optimisation of the docked pose (E_docked)
    # ------------------------------------------------------------------
    logger.info("Running constrained xTB optimisation on docked pose …")
    docked_stdout: str = run_xtb_gfn2(docked_xyz, constraint_file)
    e_docked: float = extract_energy(docked_stdout)
    logger.info("E_docked = %.12f Eh", e_docked)

    # ------------------------------------------------------------------
    # Step 4 – unconstrained optimisation of the free reference (E_free)
    # Convert the reference SDF to XYZ, preserving existing H atoms.
    # ------------------------------------------------------------------
    ref_free_xyz: str = "ref_free.xyz"
    ref_mol: Chem.Mol = next(Chem.SDMolSupplier(ref_sdf, removeHs=False))
    if ref_mol is None:
        raise IOError(f"Could not load reference SDF from {ref_sdf}")
    if ref_mol.GetNumConformers() == 0:
        raise ValueError(
            f"Reference SDF '{ref_sdf}' contains no 3-D conformer. "
            f"Generate one with RDKit EmbedMolecule or Conformator before running QUEST."
        )
    Chem.MolToXYZFile(ref_mol, ref_free_xyz)

    logger.info("Running unconstrained xTB optimisation on free reference …")
    free_stdout: str = run_xtb_gfn2(ref_free_xyz)
    e_free: float = extract_energy(free_stdout)
    logger.info("E_free  = %.12f Eh", e_free)

    # ------------------------------------------------------------------
    # Step 5 – convert to kcal mol⁻¹ and return
    # ------------------------------------------------------------------
    strain_kcal: float = (e_docked - e_free) * 627.509
    logger.info(
        "Strain energy = (%.12f - %.12f) × 627.509 = %.4f kcal mol⁻¹",
        e_docked, e_free, strain_kcal,
    )
    return strain_kcal
