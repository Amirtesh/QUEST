# QUEST: Quantum Energetics Strain Tool

QUEST is a command-line toolkit for evaluating conformational strain energies of docked ligand poses using the GFN2-xTB semiempirical quantum-mechanical method. It integrates tautomer enumeration, QM geometry optimisation, molecular docking, and strain scoring into a single end-to-end pipeline.

---

## Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Command Reference](#command-reference)
  - [quest strain](#quest-strain)
  - [quest states](#quest-states)
  - [quest dock](#quest-dock)
  - [quest run](#quest-run)
  - [quest batch](#quest-batch)
- [Output](#output)
- [Module Reference](#module-reference)
  - [Module 1: Strain Evaluator](#module-1-strain-evaluator)
  - [Module 2: State Generator](#module-2-state-generator)
  - [Module 3: Docking Engine](#module-3-docking-engine)
  - [Pipeline Orchestrator](#pipeline-orchestrator)
- [Environment Variables](#environment-variables)
- [Logging](#logging)

---

## Overview

The central scientific question QUEST addresses is: how much energetic penalty does a ligand pay by adopting the geometry required for receptor binding? A large strain energy indicates the docked conformation is far from the ligand's intrinsic energy minimum, which is an unfavourable thermodynamic signal regardless of the raw docking score.

QUEST computes this quantity as:

```
delta_E_strain = (E_docked - E_free) * 627.509   [kcal / mol]
```

where `E_docked` is the GFN2-xTB energy of the ligand locked in its docked geometry (heavy atoms frozen, hydrogens relaxed) and `E_free` is the GFN2-xTB energy of the same ligand after unconstrained optimisation from the reference geometry, both evaluated with ALPB implicit solvation (water).

---

## Architecture

QUEST is organised into three independent modules plus a master pipeline:

```
SMILES
  |
  +-- [Module 2: State Generator]
  |     Tautomer enumeration (RDKit TautomerEnumerator)
  |     3-D embedding (ETKDGv3)
  |     QM geometry optimisation (CREST iMTD-GC / GFN2)
  |     --> tautomer_N_qm.sdf
  |
  +-- [Module 3: Docking Engine]
  |     Binding-pocket detection (P2Rank)
  |     Ligand format conversion (OpenBabel SDF -> PDBQT)
  |     Molecular docking (AutoDock Vina)
  |     Pose format conversion (OpenBabel PDBQT -> SDF)
  |     --> stem_docked_poses.sdf
  |
  +-- [Module 1: Strain Evaluator]
        MCS-guided coordinate transfer (RDKit rdFMCS)
        xTB constrained optimisation (GFN2 / ALPB water)
        xTB unconstrained optimisation (GFN2 / ALPB water)
        --> strain energy in kcal / mol
```

The `quest run` and `quest batch` commands chain all three modules automatically and write a results CSV.

---

## Dependencies

### Python packages

| Package | Version | Notes |
|---------|---------|-------|
| Python | >= 3.10 | |
| click | >= 8.0 | CLI framework; installed automatically via pip |
| pandas | >= 1.5 | Results DataFrame and CSV output |
| rdkit | any recent | Cheminformatics; install via conda-forge |

### External binaries

| Binary | Purpose | Recommended installation |
|--------|---------|--------------------------|
| xtb | GFN2-xTB geometry optimisation (Modules 1 and 2) | `conda install -c conda-forge xtb` |
| crest | Conformer/minima search for QM geometries (Module 2) | `conda install -c conda-forge crest` |
| obabel | SDF/PDBQT interconversion (Module 3) | `conda install -c conda-forge openbabel` |
| vina | AutoDock Vina molecular docking (Module 3) | `conda install -c conda-forge autodock-vina` |
| prank | P2Rank binding-pocket detection (Module 3) | bundled at `tools/p2rank/prank`; see below |

### P2Rank

A copy of P2Rank is bundled at `tools/p2rank/prank`. QUEST resolves the executable in this priority order:

1. The `QUEST_P2RANK` environment variable (absolute path to `prank`).
2. The bundled copy at `<project_root>/tools/p2rank/prank`.
3. `prank` found on the system `PATH`.

If none of these is available, `quest dock` and `quest run` will raise a `FileNotFoundError` with instructions.

---

## Installation

Clone the repository, create and activate a dedicated conda environment, then install QUEST in editable mode:

```bash
git clone https://github.com/Amirtesh/QUEST.git
cd quest
conda create -n quest python=3.10
conda activate quest
conda install -c conda-forge rdkit xtb crest openbabel autodock-vina
pip install -e .
```

After installation the `quest` command is available in the activated environment.

To verify:

```bash
quest --version
```

---

## Command Reference

All commands share a global `--log-level` option and a `--version` flag:

```
quest [--log-level LEVEL] [--version] COMMAND [ARGS]...
```

`LEVEL` is one of `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` (default: `WARNING`). Pass `--log-level INFO` or use the per-command `-v` flag for detailed progress output.

---

### quest strain

Calculate the conformational strain energy of a single docked pose relative to a reference free-ligand geometry.

```
quest strain [OPTIONS] REF_SDF DOCKED_SDF
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `REF_SDF` | Reference (free) ligand SDF containing a 3-D conformer. This is the QM-minimum geometry produced by `quest states`. |
| `DOCKED_SDF` | Docked pose SDF as produced by the docking engine. |

**Options**

| Option | Default | Description |
|--------|---------|-------------|
| `-v`, `--verbose` | off | Print intermediate energies and xTB command details. Promotes log level to INFO. |
| `-o`, `--output-dir` | `.` | Directory where intermediate files (`docked_ready.xyz`, `constrain.inp`, `ref_free.xyz`) are written. |

**Output**

Prints the strain energy in kcal / mol to stdout. A positive value means the docked pose is higher in energy than the relaxed free ligand.

```
quest strain ref.sdf docked.sdf
quest strain ref.sdf docked.sdf -v --output-dir ./scratch
```

**Strain thresholds** (when `-v` is set)

| Range | Interpretation |
|-------|---------------|
| < 1.5 kcal / mol | Low strain |
| 1.5 to 5.0 kcal / mol | Moderate strain |
| > 5.0 kcal / mol | High strain |

---

### quest states

Enumerate tautomers of a ligand and find the QM-minimum geometry for each using CREST.

```
quest states [OPTIONS] SMILES
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `SMILES` | Input SMILES string. Quote it in the shell if it contains brackets, e.g. `"CC(=O)Nc1ccc(O)cc1"`. |

**Options**

| Option | Default | Description |
|--------|---------|-------------|
| `-n`, `--max-tautomers` | 3 | Maximum number of tautomers to enumerate (1 to 20). |
| `-T`, `--threads` | 4 | CPU threads passed to CREST via `-T`. |
| `-o`, `--output-dir` | `.` | Directory where `tautomer_N_qm.sdf` files are written. |
| `-v`, `--verbose` | off | Promote log level to INFO so CREST progress is visible. |

**Output**

Writes one SDF file per successfully processed tautomer, named `tautomer_0_qm.sdf`, `tautomer_1_qm.sdf`, etc. Prints the list of saved paths to stdout.

```bash
quest states "CC(=O)Nc1ccc(O)cc1"
quest states "CC(=O)Nc1ccc(O)cc1" -n 5 -T 8 -o ./states_out -v
```

---

### quest dock

Dock a single QM-optimised ligand SDF into a prepared receptor using AutoDock Vina, with optional automatic pocket detection via P2Rank.

```
quest dock [OPTIONS] QM_SDF RECEPTOR_PDBQT
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `QM_SDF` | QM-optimised ligand SDF, as produced by `quest states`. |
| `RECEPTOR_PDBQT` | Prepared receptor PDBQT file for Vina. |

**Options**

| Option | Default | Description |
|--------|---------|-------------|
| `--size SX SY SZ` | 20.0 20.0 20.0 | Docking box dimensions in Angstrom (x y z). |
| `--center CX CY CZ` | none | Docking box centre in Angstrom. If omitted, P2Rank detects the pocket automatically. |
| `--p2rank-pdb` | none | Receptor PDB for P2Rank automatic pocket detection. Required when `--center` is omitted. |
| `-e`, `--exhaustiveness` | 16 | Vina exhaustiveness (1 to 128). Higher values are more thorough but slower. |
| `-o`, `--output-dir` | `.` | Directory where the docked PDBQT and SDF are written. |
| `-v`, `--verbose` | off | Promote log level to INFO. |

**Output**

Writes `{stem}_docked.pdbqt` and `{stem}_docked_poses.sdf` to the output directory. Prints the absolute path to the SDF file.

```bash
quest dock tautomer_0_qm.sdf receptor.pdbqt --p2rank-pdb receptor.pdb
quest dock tautomer_0_qm.sdf receptor.pdbqt --center 10.5 -3.2 22.0 --size 22 22 22
quest dock tautomer_0_qm.sdf receptor.pdbqt --p2rank-pdb receptor.pdb -e 32 -v -o ./dock_out
```

---

### quest run

Run the complete QUEST pipeline for a single ligand SMILES: tautomer enumeration, QM optimisation, docking, and strain scoring.

```
quest run [OPTIONS] SMILES RECEPTOR_PDBQT
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `SMILES` | Ligand SMILES string. |
| `RECEPTOR_PDBQT` | Prepared receptor PDBQT for Vina. |

**Options**

| Option | Default | Description |
|--------|---------|-------------|
| `--p2rank-pdb` | required | Receptor PDB for P2Rank pocket detection. |
| `--size SX SY SZ` | 20.0 20.0 20.0 | Vina docking box dimensions in Angstrom. |
| `-n`, `--max-tautomers` | 3 | Maximum tautomers to enumerate. |
| `-T`, `--threads` | 4 | CPU threads for CREST. |
| `-e`, `--exhaustiveness` | 16 | Vina exhaustiveness. |
| `--strain-cutoff` | 7.0 | Strain threshold in kcal / mol for the `Viable` flag. |
| `--output-csv` | `quest_results.csv` | Path for the output results CSV. |
| `-o`, `--output-dir` | `quest_output` | Root directory for all pipeline output files. |
| `-v`, `--verbose` | off | Full pipeline progress output. |

**Pipeline stages**

1. CREST / GFN2: enumerate tautomers and QM-minimise each (Module 2).
2. P2Rank: detect the top binding pocket.
3. Vina: dock each tautomer into the pocket (Module 3).
4. xTB GFN2: calculate conformational strain per docked pose (Module 1).

**Output CSV columns**

| Column | Description |
|--------|-------------|
| `Tautomer` | Source tautomer SDF filename |
| `Pose` | Pose index within the docked SDF |
| `Vina_Affinity` | Binding affinity in kcal / mol (more negative is better) |
| `QM_Strain` | Conformational strain energy in kcal / mol |
| `Viable` | `True` if strain is below `--strain-cutoff`, else `False` |

Results are sorted by `Vina_Affinity` ascending (most negative first), then `QM_Strain` ascending.

```bash
quest run "CC(=O)Oc1ccccc1C(=O)O" receptor.pdbqt --p2rank-pdb receptor.pdb --size 20 20 20
quest run "CC(=O)Oc1ccccc1C(=O)O" receptor.pdbqt --p2rank-pdb receptor.pdb \
    --size 22 22 22 -n 5 -T 12 -e 32 --strain-cutoff 5.0 -v -o ./run1
```

---

### quest batch

Run the full QUEST pipeline over a collection of ligands from a SMILES file or a directory of SDF files.

```
quest batch [OPTIONS] INPUT_PATH RECEPTOR_PDBQT
```

**Arguments**

| Argument | Description |
|----------|-------------|
| `INPUT_PATH` | A `.smi` or `.txt` file (one SMILES per line, optionally followed by a molecule ID), or a directory of `.sdf` files. |
| `RECEPTOR_PDBQT` | Prepared receptor PDBQT for Vina. |

**Options**

| Option | Default | Description |
|--------|---------|-------------|
| `--p2rank-pdb` | required | Receptor PDB for P2Rank pocket detection. |
| `--size SX SY SZ` | 20.0 20.0 20.0 | Vina docking box dimensions in Angstrom. |
| `-n`, `--max-tautomers` | 3 | Maximum tautomers per ligand. |
| `-T`, `--threads` | 4 | CPU threads for CREST. |
| `-e`, `--exhaustiveness` | 16 | Vina exhaustiveness. |
| `--strain-cutoff` | 7.0 | Strain threshold in kcal / mol. |
| `--output-csv` | `quest_batch_results.csv` | Path for the master batch results CSV. |
| `-o`, `--output-dir` | `quest_batch_output` | Root directory. A sub-directory is created per ligand. |
| `-v`, `--verbose` | off | Full pipeline progress output. |

**SMILES file format**

```
# comment lines are ignored
CC(=O)O                    aspirin
CC(=O)Nc1ccc(O)cc1         acetaminophen
CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C
```

If no molecule ID is provided, QUEST assigns `mol_0001`, `mol_0002`, etc.

**Output**

One sub-directory per ligand under `--output-dir`, each containing the per-ligand pipeline outputs. A master CSV combining all ligands is written to `--output-csv`. The master CSV has an additional `Ligand_ID` column prepended.

```bash
quest batch ligands.smi receptor.pdbqt --p2rank-pdb receptor.pdb
quest batch ligands.smi receptor.pdbqt --p2rank-pdb receptor.pdb \
    --size 22 22 22 -n 2 -T 12 -e 32 -v -o ./batch_run
quest batch ./sdf_library/ receptor.pdbqt --p2rank-pdb receptor.pdb -v
```

---

## Output

### Intermediate files

| File | Produced by | Description |
|------|-------------|-------------|
| `tautomer_N_qm.sdf` | `quest states` | QM-minimum geometry of tautomer N |
| `p2rank-pocket/` | `quest dock` | Full P2Rank output directory with predictions CSV |
| `{stem}_docked.pdbqt` | `quest dock` | Raw Vina docked poses in PDBQT format |
| `{stem}_docked_poses.sdf` | `quest dock` | Docked poses converted to SDF by OpenBabel |
| `docked_ready.xyz` | `quest strain` | Topology-corrected XYZ of the docked pose |
| `constrain.inp` | `quest strain` | xTB heavy-atom constraint file |
| `ref_free.xyz` | `quest strain` | Reference ligand XYZ for unconstrained optimisation |

### Results CSV

The `quest run` and `quest batch` commands write a results CSV that can be opened directly in any spreadsheet application or loaded with pandas:

```python
import pandas as pd
df = pd.read_csv("quest_results.csv")
viable = df[df["Viable"]]
print(viable.sort_values("Vina_Affinity"))
```

---

## Module Reference

### Module 1: Strain Evaluator

**File:** `quest/strain_evaluator.py`

Computes the conformational strain energy of a single docked pose.

| Function | Signature | Description |
|----------|-----------|-------------|
| `transfer_topology` | `(ref_sdf, docked_sdf, output_xyz) -> int` | MCS-guided coordinate transplant from docked SDF onto the reference molecular graph; re-adds H atoms and writes an XYZ file. |
| `create_xtb_constraint` | `(filepath) -> None` | Writes an xTB constraint file that freezes heavy atoms (C, O, N, S, P, F, Cl, Br, I) with force constant 1.0 a.u. |
| `run_xtb_gfn2` | `(xyz_file, constraint_file) -> str` | Subprocess wrapper: runs `xtb --gfn 2 --alpb water --opt` in a scratch directory with a 20-minute timeout. Returns full stdout. |
| `extract_energy` | `(xtb_stdout) -> float` | Parses the `TOTAL ENERGY` line from xTB stdout and returns the value in Hartrees. |
| `calculate_strain` | `(ref_sdf, docked_sdf) -> float` | Orchestrator: runs the full strain pipeline and returns delta_E in kcal / mol. |

---

### Module 2: State Generator

**File:** `quest/state_generator.py`

Enumerates biologically relevant tautomers and obtains a QM-level minimum geometry for each.

| Function | Signature | Description |
|----------|-----------|-------------|
| `generate_tautomers` | `(smiles, max_tautomers) -> list[Mol]` | Uses RDKit `TautomerEnumerator` to enumerate and rank tautomers by stability score. |
| `embed_3d` | `(mol) -> Mol` | Adds explicit hydrogens and generates a 3-D conformer using RDKit ETKDGv3 with `randomSeed=42`. |
| `run_crest_minimum_search` | `(xyz_path, threads) -> str` | Runs CREST iMTD-GC conformer search in a temporary directory; returns the contents of `crest_best.xyz`. Environment variables `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS` are set to `1` to prevent threading conflicts. Timeout: 2 hours. |
| `update_mol_from_xyz_string` | `(mol, xyz_string) -> Mol` | Overwrites conformer coordinates from a CREST XYZ string using direct positional atom mapping. |
| `prep_ligand_states` | `(smiles, output_dir, max_tautomers, threads) -> list[str]` | Top-level orchestrator: tautomers -> embed -> CREST -> SDF. Returns a list of absolute SDF paths. |

---

### Module 3: Docking Engine

**File:** `quest/docking_engine.py`

Automates binding-pocket detection and ligand docking.

| Function | Signature | Description |
|----------|-----------|-------------|
| `_resolve_p2rank` | `(hint) -> str` | Internal: resolves the `prank` executable via env var, bundled copy, or PATH. |
| `run_p2rank` | `(pdb_path, p2rank_exec, output_dir) -> tuple[float, float, float]` | Runs P2Rank pocket prediction, parses the `_predictions.csv`, and returns `(center_x, center_y, center_z)` of the top-ranked pocket. Timeout: 10 minutes. |
| `prepare_ligand_pdbqt` | `(sdf_path, pdbqt_path) -> None` | Converts an SDF to PDBQT via OpenBabel with Gasteiger partial charges. |
| `run_vina` | `(receptor_pdbqt, ligand_pdbqt, center, size, output_pdbqt, exhaustiveness) -> str` | Runs AutoDock Vina; returns Vina stdout (contains the affinity table). Timeout: 1 hour. |
| `convert_vina_output_to_sdf` | `(vina_out_pdbqt, output_sdf) -> None` | Converts a multi-pose Vina PDBQT to SDF via OpenBabel. |
| `dock_ligand` | `(qm_sdf, receptor_pdbqt, size, center, p2rank_pdb, ...) -> str` | Orchestrator: pocket detection -> PDBQT conversion -> Vina docking -> SDF conversion. Returns the absolute path to the docked SDF. |

---

### Pipeline Orchestrator

**File:** `quest/pipeline.py`

| Function | Signature | Description |
|----------|-----------|-------------|
| `_extract_vina_affinity` | `(pose) -> float or None` | Extracts the Vina binding affinity from an RDKit Mol read from a Vina-output SDF, trying `minimizedAffinity`, `REMARK`, and `_Name` properties in order. |
| `run_quest_pipeline` | `(smiles, receptor_pdb, receptor_pdbqt, size, ...) -> DataFrame` | Runs the four-stage pipeline for a single ligand and returns a results DataFrame. Sorts by Vina affinity then strain energy. |
| `run_quest_batch` | `(input_path, receptor_pdb, receptor_pdbqt, size, ...) -> DataFrame` | Iterates over a SMILES file or SDF directory, calls `run_quest_pipeline` per molecule, and concatenates results into a master DataFrame with a `Ligand_ID` column. |

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `QUEST_P2RANK` | Absolute path to the `prank` executable. Overrides the bundled copy and PATH lookup. |

---

## Logging

QUEST uses Python's standard `logging` module. The root logger is configured by the `--log-level` option on the `quest` group or by the per-command `-v` / `--verbose` flag (which promotes the level to `INFO`).

To capture full debug output including all subprocess commands and intermediate energies:

```bash
quest --log-level DEBUG run "CC(=O)O" receptor.pdbqt --p2rank-pdb receptor.pdb
```

Log messages follow the format:

```
HH:MM:SS  <logger-name>                  <LEVEL>   <message>
```
