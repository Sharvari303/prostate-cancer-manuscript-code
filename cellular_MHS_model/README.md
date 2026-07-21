# Prostate Cancer Hybrid (Multiscale) Model

A hybrid **ODE + Boolean** model of prostate cancer cell-fate decisions. Continuous
kinetic pathways (ErbB-mediated Ras–MAPK / PI3K–AKT and Androgen Receptor signaling)
are simulated in **COPASI**, while an upstream p53/cell-cycle regulatory network is
simulated as a **discrete Boolean network**. The two layers are coupled: at every
Boolean step the ODE steady states are discretized and injected into the Boolean
network, and the resulting Boolean node states re-parameterize the next ODE run.
The end product is a set of **cell-fate probabilities** (death / proliferation /
senescence) aggregated over many Boolean initial conditions.

> Model by Dr. Ravi Radhakrishnan and Alok Ghosh.

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Environment Setup](#environment-setup)
3. [Quick Start](#quick-start)
4. [Running on an HPC Cluster (Slurm)](#running-on-an-hpc-cluster-slurm)
5. [Command-Line Options](#command-line-options)
6. [How the Simulation Works](#how-the-simulation-works)
7. [Input Files](#input-files)
8. [Output: Where Everything Gets Stored](#output-where-everything-gets-stored)
9. [Reproducing the Manuscript Figures (Parameter Sweeps)](#reproducing-the-manuscript-figures-parameter-sweeps)
10. [File & Directory Reference](#file--directory-reference)

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| **Python 2.7** (≥ 2.6) |  |
| **COPASI (`CopasiSE`)** | Bundled — no separate install needed. See [Environment Setup](#environment-setup). |
| **NumPy** |   |
| **Matplotlib** (≥ 1.5) |  | 

The COPASI command-line engine `CopasiSE` is provided in-repo for both architectures:

```
Copasi_Executables/Linux_64Bit/CopasiSE      <- used on 64-bit systems (default)
Copasi_Executables/Linux_32Bit/CopasiSE      <- used on 32-bit systems
```

The main script invokes `CopasiSE` by bare name, so **that directory must be on your
`PATH`** at run time (the Slurm script and the venv setup both handle this for you).

---

## Environment Setup

On the **HPC cluster**, create the virtual environment once (module names may differ
by site — adjust to whatever provides Python 2.7):

```bash
module load python/2.7.18
python2 -m virtualenv venv_pca
source venv_pca/bin/activate
pip install -r requirements.txt
```

Then, in every future session:

```bash
module load python/2.7.18
source venv_pca/bin/activate
export PATH=$PWD/Copasi_Executables/Linux_64Bit:$PATH   # make CopasiSE findable
```

Verify the engine is reachable:

```bash
which CopasiSE      # should print a path inside Copasi_Executables/Linux_64Bit
```

---

## Quick Start

From inside the model directory, with the environment activated and `CopasiSE` on `PATH`:

```bash
python combined_ode_boolean.py -p Model_Input.json &> Output.log
```

- `-p` runs the different Boolean initial conditions **in parallel** (8 cores by default).
- `Model_Input.json` is the master input file (all ODE + Boolean model definitions).
- All results are written to a new `Output_Data/` directory created in the working directory.
- When it finishes, the script prints the cell-fate summary to stdout, e.g.:

```
Simulation Results:
cell_death,cell_growth,cell_senescence
0.458333,0.166667,0.375000
```

A single run with the default `Model_Input.json` (8 initial conditions × 2 Boolean
steps) takes on the order of ~7–8 minutes **per initial condition**; with `-p` on 8
cores those run concurrently.

---

## Running on an HPC Cluster (Slurm) 

Need to adapt for use. The provided `Job_Submit_Slurm.sh` bundles everything (module load, venv activation,
`PATH` export, and the run). Edit the `#SBATCH` directives (partition/queue, account,
cores, wall time) and the `module load` line to match your cluster before submitting:

```bash
sbatch Job_Submit_Slurm.sh
```

---

## Command-Line Options

Run `python combined_ode_boolean.py --help` to see all options:

| Flag | Long form | Default | Purpose |
|------|-----------|---------|---------|
| `-o DIR` | `--output-dirname` | `Output_Data` | Output directory (absolute, or relative to the model folder). |
| `-p` | `--parallel-run` | off | Run the initial conditions in parallel across cores. |
| `-n N` | `--no-of-cores` | 8 (when `-p`) | Number of worker cores for the parallel pool. |
| `-r` | `--reinitialize-ode` | off | **"No-memory" mode** — reset the ODE modules to base state at every Boolean step instead of carrying concentrations forward. Used to represent cell-cycle variation. |
| `-M FILE` | `--mirtar-file` | `hsa_MTI.csv` | miRTarBase CSV for patient-specific miRNA effects (not used in the base study). |

The positional argument (`Model_Input.json` above) is the required input JSON.

---

## How the Simulation Works

The pipeline is orchestrated by `combined_ode_boolean.py`, which drives these modules:

```
combined_ode_boolean.py   <- entry point, CLI, parallel dispatch, post-processing
  └─ mwe_base.Combined_Run       <- coupled ODE<->Boolean loop for ONE initial condition
       ├─ CopasiRun.OdeModule    <- runs one COPASI module, discretizes to node states
       │    └─ interface.CopasiState  <- reads/edits/writes .cps (COPASI XML) files
       │    └─ CopasiSE           <- external COPASI engine (subprocess)
       └─ Truth_Table.Calculate_State  <- one synchronous Boolean update step
  └─ mwe_collector.Collect_Data  <- aggregates all conditions into cell-fate stats
```

Step by step:

1. **Load & prep.** `Model_Input.json` is read, `Output_Data/` is created, and the base
   COPASI file is copied to `ErbB_Ras-MAPK_PI3K-AKT_adjMirna.cps` inside the output dir.
   The list of Boolean **initial conditions** to sweep comes from `init_iterator` in the
   JSON (8 states in the default file, each a bit-vector over the *Unconstrained_Nodes*).

2. **Dispatch.** With `-p`, the initial conditions are distributed to a
   `multiprocessing.Pool` (8 workers by default). Without `-p`, they run serially.

3. **Per-condition coupled loop** (`Combined_Run`), repeated for `STEPS` iterations:
   - Run the **MAPK/PI3K COPASI module**, then the **AR COPASI module**, via `OdeModule`.
     Each call edits a temporary `.cps`, shells out to `CopasiSE`, and reads back the
     time-course CSV.
   - **Discretize** each interface species against `thresh` in the JSON: below 0.33×peak
     → node `0`, above 0.6×peak → `1`, in between → `2` (ambiguous, left unchanged).
   - Feed those discrete states into the **Boolean network** (`Calculate_State`) for one
     synchronous update, using the weighted interaction rules in `p53_network`.
   - Feed the new Boolean states back to re-parameterize the next ODE run (e.g. TP53, AR,
     ACT initial concentrations). With `-r`, the ODE modules are reset each step instead.
   - Accumulate node trajectories and state-transition probabilities.

4. **Post-processing** (`Collect_Data`). After all conditions finish, every
   `Output_Data*.pickle` is re-read and reduced to **cell-fate counts and fractions**:
   - **death** — CASP9 active while CDK2 inactive;
   - **proliferation** — high PSA (> 50) or CDK2 active without CASP9;
   - **senescence** — everything else.

---

## Input Files

| File | Role |
|------|------|
| **`Model_Input.json`** | Master input. Contains the Boolean network (`p53_network`), the ODE↔Boolean coupling config, thresholds, initial-condition sweep (`init_iterator`), number of Boolean steps (`STEPS`), COPASI file mapping, chemo-drug settings, and file paths. This is the file you edit to change a run. |
| **`Input_File_Template.json`** | A stripped-down template/example of the input format. |
| **`ErbB_Ras-MAPK_PI3K-AKT.cps`** | COPASI model of the ErbB-mediated Ras–MAPK and PI3K–AKT pathway ("mapk" module). |
| **`AR_Base_Interfaces.cps`** | COPASI model of the Androgen Receptor signaling pathway ("ar" module). |

### Key fields inside `Model_Input.json`

- **`p53_network`** — the Boolean network. Each node has a `base`, an `initial_state`,
  and weighted `input_nodes` (`[node, weight]` pairs; sign = activation/inhibition).
- **`init_iterator`** — the list of Boolean initial conditions to sweep; **each entry is
  one independent simulation.** More entries ⇒ longer run, finer cell-fate statistics.
- **`STEPS`** — number of coupled ODE↔Boolean iterations per initial condition.
- **`Constrained_Nodes` / `Unconstrained_Nodes`** — which nodes are driven by the ODE
  layer vs. which are varied by `init_iterator`. This is the place where DEG relevant
  information for patient cohorts from TCGA or EUREKA1 patient data are used.
- **`thresh`** — concentration thresholds used to discretize ODE species into Boolean
  node states.
- **`changes`** — per-module COPASI edits applied before each run (e.g. simulation
  `Duration`, `StepSize`, and select initial species values such as PTEN or Testosterone).
- **`Copasifiles`** — maps each module to its **Original** `.cps`, a per-run **Temporary**
  `.cps`, the **Output** CSV name, and the ODE steady-state JSON name.

---

## Output: Where Everything Gets Stored

Everything goes into the output directory — **`Output_Data/` by default**, or whatever
you pass with `-o`. Files whose names contain a bit-string (e.g. `...000100111001100`)
are **per-initial-condition**; that suffix is the Boolean initial-condition fingerprint.

### Result files 

| File | Format | Contents |
|------|--------|----------|
| **`cell_fates.csv`** | CSV | The headline result: one row of `cell_death,cell_growth,cell_senescence` fractions. |
| **`Cell_Fate.json`** | JSON | Same result with raw counts **and** fractions, e.g. `{"cell_death":[11,0.458], "cell_proliferation":[4,0.167], "cell_senescence":[9,0.375]}`. |
| **`ckr_total`** | text | The cell-kill fraction as a single number (convenience for sweeps). |


---

## Reproducing the Manuscript Figures (Parameter Sweeps)

The instructions above run a **single instance** of the model. The manuscript figures
come from a **sweep**: many independent copies of this base model, each with a different

- **Testosterone level** (edit the `medium@Testosterone_Buffered` value under
  `changes.ar` in `Model_Input.json`),
- **Growth-factor level** (edit the corresponding initial species in `changes`),
- **Memory condition** (add/remove the `-r` flag).

Each copy is run independently on the cluster (one Slurm job each). After all jobs
finish, collect the `cell_fates.csv` / `Cell_Fate.json` values from every run and plot
them. The preprocessing scripts that generated the PBS/Slurm sweep directories are
site-specific and are **not** included here (see the legacy `Job_Submit.sh` for the
directory-per-case pattern that was used).

The `collected_csvs/` directory holds spreadsheets of previously collected sweep results
(e.g. net growth rates, PTEN expressions, patient data) used downstream for figures.

---

## File & Directory Reference

| Path | What it is |
|------|-----------|
| `combined_ode_boolean.py` | **Main entry point.** CLI, parallel dispatch, post-processing. |
| `mwe_base.py` (`Combined_Run`) | Coupled ODE↔Boolean loop for one initial condition. |
| `CopasiRun.py` (`OdeModule`) | Runs one COPASI module and discretizes its output. |
| `mwe_collector.py` (`Collect_Data`) | Aggregates all conditions into cell-fate stats. |
| `interface.py` (`CopasiState`) | Reads/edits/writes COPASI `.cps` (XML) files. |
| `Truth_Table.py` (`Calculate_State`) | One synchronous Boolean-network update. |
| `interpret_link.py` | miRNA target interpretation (miRTarBase). |
| `tools.py`, `ordereddict.py` | Helpers|
| `Model_Input.json` | Master input file. |
| `Input_File_Template.json` | Input format template. |
| `ErbB_Ras-MAPK_PI3K-AKT.cps` | MAPK/PI3K COPASI model. |
| `AR_Base_Interfaces.cps` | AR COPASI model. |
| `Copasi_Executables/` | Bundled `CopasiSE` engines. |
| `Job_Submit_Slurm.sh` | Slurm submission script (adapt `#SBATCH` directives to your cluster). |
| `Job_Submit.sh` | Legacy PBS submission script. |
| `requirements.txt` | Python 2 dependencies. |
| `venv_pca/` | Python 2.7 virtual environment. |
| `Output_Data/` | Default output directory (see [Output](#output-where-everything-gets-stored)). |
| `collected_csvs/` | Collected sweep results / downstream data. |
| `run.log` | Example console log from a previous run. |
```
