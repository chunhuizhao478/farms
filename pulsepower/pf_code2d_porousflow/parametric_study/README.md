# Porousflow Parametric Study Scripts Documentation

This directory contains automated scripts for running parametric studies on porousflow simulations.

## Table of Contents

- [Overview](#overview)
- [Scripts](#scripts)
  - [1. generate_parametric_cases.py](#1-generate_parametric_casespy)
  - [2. run_static_solve_and_update.py](#2-run_static_solve_and_updatepy)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)

---

## Overview

These scripts automate the workflow for parametric studies involving:
1. **Case generation** - Creating multiple simulation cases with varying parameters
2. **Static solve execution** - Running steady-state simulations
3. **Energy extraction** - Extracting energy values from simulation results
4. **File updates** - Automatically updating dynamic simulation files with static results

### Workflow

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Generate parametric cases                                │
│    → Creates folders with modified input files              │
└─────────────────────────┬───────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. Run static solve                                         │
│    → Executes static_solve.i for each case                  │
│    → Produces static_solve_csv.csv output                   │
└─────────────────────────┬───────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. Extract energy values and update elasticity.i            │
│    → Reads CSV output                                       │
│    → Updates elasticity.i with static energy values         │
└─────────────────────────────────────────────────────────────┘
                          ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. Ready for dynamic simulation                             │
│    → Run elasticity.i (with correct initial conditions)     │
└─────────────────────────────────────────────────────────────┘
```

---

## Scripts

### 1. `generate_parametric_cases.py`

**Purpose:** Generates parametric study cases by creating folders and modified input files.

#### Features

- Creates two types of parametric studies:
  - **Study 1:** Confinement pressure variations (1, 5, 10 MPa)
  - **Study 2:** Initial pore pressure variations (0.0965, 2, 4 MPa) with fixed 5 MPa confinement
- Supports multiple domain sizes (1x, 2x, 5x)
- Modifies three input files per case:
  - `static_solve.i` - Static equilibrium simulation
  - `elasticity.i` - Dynamic simulation
  - `fracture.i` - Fracture mechanics sub-app

#### Usage

```bash
# Navigate to the parametric study directory
cd /path/to/parametric_study

# Generate all parametric cases
python3 generate_parametric_cases.py
```

#### Output

The script generates **18 case folders** in total:

**Study 1 - Confinement Variations (9 cases):**
```
case_cf1_domain1x/
case_cf1_domain2x/
case_cf1_domain5x/
case_cf5_domain1x/
case_cf5_domain2x/
case_cf5_domain5x/
case_cf10_domain1x/
case_cf10_domain2x/
case_cf10_domain5x/
```

**Study 2 - Pore Pressure Variations (9 cases):**
```
case_pp0d0965_domain1x/
case_pp0d0965_domain2x/
case_pp0d0965_domain5x/
case_pp2_domain1x/
case_pp2_domain2x/
case_pp2_domain5x/
case_pp4_domain1x/
case_pp4_domain2x/
case_pp4_domain5x/
```

#### Configuration

Edit the **CONFIGURATION** section at the top of the script:

```python
# Confinement pressures to study
CONFINEMENT_PRESSURES = {
    'cf1': 1e6,   # 1 MPa
    'cf5': 5e6,   # 5 MPa
    'cf10': 10e6  # 10 MPa
}

# Initial pore pressures to study
INITIAL_PORE_PRESSURES = {
    'pp0d0965': 0.0965e6,  # 0.0965 MPa
    'pp2': 2e6,            # 2 MPa
    'pp4': 4e6             # 4 MPa
}

# Fixed confinement for pore pressure study
FIXED_CONFINEMENT_FOR_PP_STUDY = 5e6  # 5 MPa

# Domain sizes
DOMAIN_SIZES = {
    'domain1x': {
        'mesh_file': '../../../2dmeshfile/fieldscale_test1_2d.msh',
        'coord': '0.01 0.01 0'
    },
    'domain2x': {
        'mesh_file': '../../../2dmeshfile/fieldscale_test1_2d_extend2x.msh',
        'coord': '0.02 0.02 0'
    },
    'domain5x': {
        'mesh_file': '../../../2dmeshfile/fieldscale_test1_2d_extend5x.msh',
        'coord': '0.06 0.06 0'
    }
}
```

#### Parameters Modified

For each case, the script modifies:

| File | Parameters Modified |
|------|---------------------|
| `static_solve.i` | • `confinement_pressure`<br>• `initial_pore_pressure`<br>• Mesh file path<br>• Fixed point coordinate |
| `elasticity.i` | • `confinement_pressure`<br>• Mesh file path<br>• Fixed point coordinate |
| `fracture.i` | • Mesh file path<br>• Fixed point coordinate |

---

### 2. `run_static_solve_and_update.py`

**Purpose:** Runs static solve simulations and automatically updates elasticity.i files with extracted energy values.

#### Features

- Executes `static_solve.i` using MPI
- Extracts energy values from CSV output
- Updates first 3 lines of `elasticity.i` automatically
- Supports batch processing of multiple cases
- Flexible case selection (single, pattern, all)
- Performance-optimized (no output buffering)

#### Usage

##### Basic Usage

```bash
# Process a specific case
python3 run_static_solve_and_update.py --case case_cf1_domain1x

# Process all confinement cases
python3 run_static_solve_and_update.py --pattern "case_cf*"

# Process all pore pressure cases
python3 run_static_solve_and_update.py --pattern "case_pp*"

# Process ALL cases
python3 run_static_solve_and_update.py --all
```

##### Advanced Options

```bash
# Dry run - extract values without updating files
python3 run_static_solve_and_update.py --case case_cf1_domain1x --dry-run

# Skip simulation - use existing CSV files
python3 run_static_solve_and_update.py --case case_cf1_domain1x --skip-simulation

# Quiet mode - suppress MOOSE output for faster execution
python3 run_static_solve_and_update.py --case case_cf1_domain1x --quiet

# Specify number of MPI processes
python3 run_static_solve_and_update.py --case case_cf1_domain1x --np 16

# Combine options
python3 run_static_solve_and_update.py --all --quiet --np 16
```

#### Command-Line Options

| Option | Short | Description |
|--------|-------|-------------|
| `--case CASE` | | Process a specific case folder |
| `--pattern PATTERN` | | Process all cases matching glob pattern |
| `--all` | | Process all case folders |
| `--dry-run` | | Extract values but don't update elasticity.i |
| `--skip-simulation` | | Skip running simulation (use existing CSV) |
| `--quiet` | `-q` | Suppress MOOSE output for speed |
| `--np NP` | | Number of MPI processes (default: 8) |
| `--help` | `-h` | Show help message |

#### What It Does

1. **Runs simulation:**
   ```bash
   mpirun -np 8 ./farms-opt -i pulsepower/pf_code2d_porousflow/parametric_study/case_cf1_domain1x/static_solve.i
   ```

2. **Reads CSV output:**
   ```
   case_cf1_domain1x/static_solve_csv.csv
   ```

   Extracts from last row (steady-state values):
   - `fluid_elastic_energy_total_static`
   - `solid_elastic_energy_static` (maps to `solid_elastic_energy_total_static`)
   - `full_input_energy_static`

3. **Updates elasticity.i:**
   ```
   fluid_elastic_energy_total_static = 8.081664e-05
   solid_elastic_energy_total_static = 5.408951e-03
   full_input_energy_static = 5.489768e-03
   ```

#### Configuration

Edit the **CONFIGURATION** section:

```python
# Number of MPI processes to use
MPI_PROCESSES = 8

# Executable name (relative to project root)
EXECUTABLE = "./farms-opt"

# CSV output filename
CSV_FILENAME = "static_solve_csv.csv"

# Column name mapping (CSV → elasticity.i)
CSV_COLUMN_NAMES = {
    'fluid_elastic_energy_total_static': 'fluid_elastic_energy_total_static',
    'solid_elastic_energy_total_static': 'solid_elastic_energy_static',
    'full_input_energy_static': 'full_input_energy_static'
}
```

#### Output Example

```
================================================================================
Static Solve Runner and Elasticity File Updater
================================================================================

Script directory: /path/to/parametric_study
Project root: /path/to/farms_cdms
Executable: ./farms-opt
MPI processes: 8

Found 1 case(s) to process:
  - case_cf1_domain1x

================================================================================
Starting processing...
================================================================================

================================================================================
Processing: case_cf1_domain1x
================================================================================

  Running: mpirun -np 8 ./farms-opt -i pulsepower/.../static_solve.i
  Working directory: /path/to/farms_cdms
  MOOSE output will be displayed below...
  ----------------------------------------------------------------------------
  [MOOSE output appears here...]
  ----------------------------------------------------------------------------
  ✓ Static solve completed successfully

  Extracting values from CSV...
  ✓ Extracted values:
    fluid_elastic_energy_total_static = 8.081664e-05
    solid_elastic_energy_total_static = 5.408951e-03
    full_input_energy_static = 5.489768e-03

  Updating elasticity.i...
  ✓ Updated elasticity.i with new values:
    fluid_elastic_energy_total_static = 8.081664e-05
    solid_elastic_energy_total_static = 5.408951e-03
    full_input_energy_static = 5.489768e-03

================================================================================
SUMMARY
================================================================================

Total cases: 1
Successful: 1
Failed: 0

✓ Successful cases:
  - case_cf1_domain1x

================================================================================
```

---

## Quick Start

### Complete Workflow Example

```bash
# 1. Generate all parametric cases
python3 generate_parametric_cases.py

# 2. Process all cases (run static solve and update elasticity.i)
python3 run_static_solve_and_update.py --all --quiet --np 16

# 3. Run dynamic simulations (example for one case)
cd case_cf1_domain1x
mpirun -np 16 ../../farms-opt -i elasticity.i
```

### Process Specific Study

```bash
# Generate cases
python3 generate_parametric_cases.py

# Process only confinement pressure study cases
python3 run_static_solve_and_update.py --pattern "case_cf*" --quiet

# Process only pore pressure study cases
python3 run_static_solve_and_update.py --pattern "case_pp*" --quiet
```

### Test Before Running

```bash
# Generate cases
python3 generate_parametric_cases.py

# Test one case with dry-run
python3 run_static_solve_and_update.py --case case_cf1_domain1x --dry-run

# If successful, process all
python3 run_static_solve_and_update.py --all
```

---

## Configuration

### Directory Structure

```
parametric_study/
├── README.md                          # This file
├── generate_parametric_cases.py       # Case generation script
├── run_static_solve_and_update.py    # Static solve runner
├── original_files/                    # Template files (DO NOT DELETE)
│   └── case_name/
│       ├── elasticity.i
│       ├── fracture.i
│       └── static_solve.i
├── case_cf1_domain1x/                 # Generated case folder
│   ├── elasticity.i
│   ├── fracture.i
│   ├── static_solve.i
│   └── static_solve_csv.csv          # Output (after running)
├── case_cf1_domain2x/
├── case_cf1_domain5x/
└── ... (more case folders)
```

### Requirements

- Python 3.6+
- MOOSE framework (`farms-opt` executable)
- MPI (for parallel execution)
- CSV module (Python standard library)

### Template Files

The scripts reference template files from:
```
original_files/case_name/
```

**⚠️ DO NOT DELETE THIS DIRECTORY** - It contains the base templates for all case generation.

---

## Troubleshooting

### Common Issues

#### 1. CSV file not found

**Error:**
```
✗ ERROR: CSV file not found: .../static_solve_csv.csv
```

**Solution:**
- Verify static solve completed successfully
- Check CSV output filename in MOOSE input file
- Ensure output wasn't redirected to different location

#### 2. Executable not found

**Error:**
```
✗ ERROR: Executable not found: .../farms-opt
```

**Solution:**
- Build MOOSE application: `make -j8`
- Verify executable location in script configuration
- Update `EXECUTABLE` path in script

#### 3. Column names don't match

**Error:**
```
✗ WARNING: Column 'solid_elastic_energy_static' not found in CSV
Available columns: ...
```

**Solution:**
- Check postprocessor names in `static_solve.i`
- Update `CSV_COLUMN_NAMES` mapping in script
- Verify CSV header matches expected names

#### 4. Slow execution through Python

**Problem:** Script runs slower than direct mpirun

**Solution:**
- Use `--quiet` flag to suppress output buffering
- Increase MPI processes: `--np 16`
- Example: `python3 run_static_solve_and_update.py --all --quiet --np 16`

#### 5. Permission denied

**Error:**
```
Permission denied: .../generate_parametric_cases.py
```

**Solution:**
```bash
chmod +x generate_parametric_cases.py
chmod +x run_static_solve_and_update.py
```

### Performance Tips

1. **Use quiet mode for batch processing:**
   ```bash
   python3 run_static_solve_and_update.py --all --quiet
   ```

2. **Increase MPI processes for faster simulation:**
   ```bash
   python3 run_static_solve_and_update.py --all --np 16 --quiet
   ```

3. **Process cases in parallel manually:**
   ```bash
   # Terminal 1
   python3 run_static_solve_and_update.py --pattern "case_cf1*" --quiet &

   # Terminal 2
   python3 run_static_solve_and_update.py --pattern "case_cf5*" --quiet &

   # Terminal 3
   python3 run_static_solve_and_update.py --pattern "case_cf10*" --quiet &
   ```

4. **Skip already-run simulations:**
   ```bash
   # If CSV already exists, just update elasticity.i
   python3 run_static_solve_and_update.py --all --skip-simulation
   ```

### Debugging

#### Enable verbose output

```bash
# See all MOOSE output (default)
python3 run_static_solve_and_update.py --case case_cf1_domain1x

# Dry run to test extraction without running simulation
python3 run_static_solve_and_update.py --case case_cf1_domain1x --dry-run --skip-simulation
```

#### Check generated files

```bash
# After generating cases, verify files were created
ls -la case_cf1_domain1x/

# Check first few lines of elasticity.i
head -n 10 case_cf1_domain1x/elasticity.i

# Compare with template
diff original_files/case_name/elasticity.i case_cf1_domain1x/elasticity.i
```

#### Verify CSV content

```bash
# Check CSV output from static solve
cat case_cf1_domain1x/static_solve_csv.csv

# Should show:
# time,fluid_elastic_energy_total_static,full_energy_static,full_input_energy_static,solid_elastic_energy_static
# 0,0,0,-0,0
# 1,8.0816643660968e-05,0.005489767560866,0.0054897675608653,0.0054089509172051
```

---

## Best Practices

1. **Always keep template files safe:**
   - Never modify files in `original_files/case_name/`
   - Make a backup before major changes

2. **Test with one case first:**
   ```bash
   python3 run_static_solve_and_update.py --case case_cf1_domain1x --dry-run
   ```

3. **Use version control:**
   ```bash
   git add generate_parametric_cases.py run_static_solve_and_update.py
   git commit -m "Add parametric study scripts"
   ```

4. **Document parameter changes:**
   - Keep a log of which parameters were varied
   - Note any manual modifications to generated cases

5. **Check results systematically:**
   ```bash
   # Verify all cases completed
   python3 run_static_solve_and_update.py --all --skip-simulation | grep "Total cases"
   ```

---

## Additional Resources

- **MOOSE Documentation:** https://mooseframework.inl.gov/
- **MPI Tutorial:** https://mpitutorial.com/
- **Python subprocess:** https://docs.python.org/3/library/subprocess.html

---

## Author & License

**Author:** Auto-generated
**Date:** 2025-11-08
**Version:** 1.0

For questions or issues, please contact the repository maintainer.

---

**Last Updated:** 2025-11-08
