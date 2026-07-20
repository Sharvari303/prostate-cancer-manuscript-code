# Prostate Cancer Agent-Based Model (PhysiCell)

Agent-based model of prostate tumor evolution under androgen deprivation therapy (ADT), simulating competition between PTEN-normal (androgen-sensitive, S) and PTEN-deleted (androgen-resistant, R) cells. Built using PhysiCell v1.14.0. (https://physicell.org/)

---

## Cell Types

| Cell Type | Description |
|---|---|
| PTEN_normal (S) | Androgen-sensitive cells; growth depends on testosterone via Michaelis-Menten kinetics |
| PTEN_deleted (R) | Androgen-resistant cells; reduced testosterone dependence |

Both cell types share identical mechanical parameters (adhesion, repulsion, motility) — differentiation arises through cohort-specific growth rate parameters.

---

## Microenvironment Substrates

| Substrate | Diffusion coefficient | Decay rate | Boundary condition |
|---|---|---|---|
| Oxygen | 100 µm²/min | 0 1/min | 38 mmHg (Dirichlet, all faces) |
| Testosterone | 1 µm²/min | 0 1/min | 8 ng/ml (Dirichlet, all faces) |

---

## Clinical Cohorts

Simulations are parameterized by three clinical cohorts, each with fitted growth rate parameters for S and R cells:

| Cohort | Description | cohort value in XML |
|---|---|---|
| BR | Biochemical Recurrent | 0 |
| TR | Imaging-based Recurrent | 1 |
| CTR | Control | 2 |

Growth rates, apoptosis rates, and testosterone sensitivity (Michaelis-Menten parameters m, n, p) differ across cohorts. All cohort parameter values are defined in `config/PhysiCell_settings.xml`. The masterlist CSV specifies which cohort and run configuration to use, and `ABMruns_updatexml.py` updates the XML accordingly before each run.

---

## Two Simulation Conditions

`ABM_unconstrained/` allows cells to grow freely up to the full ±250 µm mesh boundary. `ABM_densepacking/` restricts cell movement to ±125 µm, mimicking growth within a spatially confined compartment. All other parameters are identical between the two conditions.

---

## How to Reproduce a Simulation

### 1. Get PhysiCell v1.14.0
```bash
git clone https://github.com/MathCancer/PhysiCell.git
cd PhysiCell
git checkout v1.14.0
```

### 2. Copy files from this repo into your PhysiCell directory
For unconstrained runs:
```bash
cp ABM_unconstrained/config/PhysiCell_settings.xml PhysiCell/config/
cp ABM_unconstrained/config/cell_rules.csv PhysiCell/config/
cp ABM_unconstrained/custom_modules/custom.cpp PhysiCell/custom_modules/
cp ABM_unconstrained/custom_modules/custom.h PhysiCell/custom_modules/
cp ABM_unconstrained/core/PhysiCell_standard_models.cpp PhysiCell/core/
```
For spatially confined runs, use `ABM_densepacking/` instead.

### 3. Compile
```bash
cd PhysiCell
make
```

### 4. Configure a run from the masterlist
```bash
python3 ABMruns_updatexml.py /path/to/PhysiCell/config/PhysiCell_settings.xml ABMruns_masterlist_prostatecancer.csv RUN_ID
```
Replace `RUN_ID` with any integer from `ABMruns_masterlist_prostatecancer.csv`. Both files are in `runs/`.

### 5. Run
```bash
cd PhysiCell
./project
```
Or submit via SLURM using the scripts in `runs/`.

---

## Repository Structure

```
ABM_unconstrained/                        # Unconstrained simulation files
  config/
    PhysiCell_settings.xml
    cell_rules.csv
  custom_modules/
    custom.cpp
    custom.h
  core/
    PhysiCell_standard_models.cpp

ABM_densepacking/                         # Spatially confined simulation files
  config/
    PhysiCell_settings.xml
    cell_rules.csv
  custom_modules/
    custom.cpp
    custom.h
  core/
    PhysiCell_standard_models.cpp         # double boundary = 125

runs/
  ABMruns_masterlist_prostatecancer.csv   # Full parameter table for crowding runs
  ABMruns_masterlist_prostatecancer_v2.csv   # Full parameter table for adhesion-motility, androgen-uptake runs
  ABMruns_updatexml.py                    # Updates PhysiCell XML config for a given Run_ID
  README.md #csv descriptions

analysis/
  ABMruns_PCa_dataanalysis.py             # Main data analysis script

README.md
```
