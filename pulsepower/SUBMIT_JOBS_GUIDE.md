# Job Submission Guide for Parametric Studies

This guide explains how to submit SLURM jobs for both pure solid and porous flow parametric studies.

## Directory Structure

```
pulsepower/
├── pf_code2d_puresolid/
│   └── parametric_study/
│       ├── submit_all_jobs.py      # Submits ONLY pure solid cases
│       ├── submit_all_jobs.sh      # Submits ONLY pure solid cases
│       ├── README_SUBMIT_JOBS.md   # Documentation
│       └── case_*/                 # 9 pure solid cases (cf1, cf5, cf10)
│
└── pf_code2d_porousflow/
    └── parametric_study/
        ├── submit_all_jobs.py      # Submits ONLY porous flow cases
        ├── submit_all_jobs.sh      # Submits ONLY porous flow cases
        ├── README_SUBMIT_JOBS.md   # Documentation
        └── case_*/                 # 18 porous flow cases (cf* + pp*)
```

## Important: Script Isolation

**Each script ONLY operates on files in its own directory:**

- Scripts in `pf_code2d_puresolid/parametric_study/` will **ONLY** find and submit pure solid cases
- Scripts in `pf_code2d_porousflow/parametric_study/` will **ONLY** find and submit porous flow cases

This is automatic - the scripts use their own location to determine which cases to process.

## Case Overview

### Pure Solid Cases (9 cases)
Located in: `pf_code2d_puresolid/parametric_study/`

- `case_cf1_domain1x`, `case_cf1_domain2x`, `case_cf1_domain5x`
- `case_cf5_domain1x`, `case_cf5_domain2x`, `case_cf5_domain5x`
- `case_cf10_domain1x`, `case_cf10_domain2x`, `case_cf10_domain5x`

### Porous Flow Cases (18 cases)
Located in: `pf_code2d_porousflow/parametric_study/`

**Compression factor cases:**
- `case_cf1_domain1x`, `case_cf1_domain2x`, `case_cf1_domain5x`
- `case_cf5_domain1x`, `case_cf5_domain2x`, `case_cf5_domain5x`
- `case_cf10_domain1x`, `case_cf10_domain2x`, `case_cf10_domain5x`

**Pore pressure cases:**
- `case_pp0d0965_domain1x`, `case_pp0d0965_domain2x`, `case_pp0d0965_domain5x`
- `case_pp2_domain1x`, `case_pp2_domain2x`, `case_pp2_domain5x`
- `case_pp4_domain1x`, `case_pp4_domain2x`, `case_pp4_domain5x`

## Usage

### For Pure Solid Cases

```bash
# Navigate to pure solid directory
cd /path/to/farms_cdms/pulsepower/pf_code2d_puresolid/parametric_study

# Validate all pure solid cases
./submit_all_jobs.py --validate

# Submit all pure solid cases
./submit_all_jobs.py --yes

# Submit only specific compression factor
./submit_all_jobs.py --pattern "case_cf5*" --yes
```

### For Porous Flow Cases

```bash
# Navigate to porous flow directory
cd /path/to/farms_cdms/pulsepower/pf_code2d_porousflow/parametric_study

# Validate all porous flow cases
./submit_all_jobs.py --validate

# Submit all porous flow cases
./submit_all_jobs.py --yes

# Submit only compression factor cases
./submit_all_jobs.py --pattern "case_cf*" --yes

# Submit only pore pressure cases
./submit_all_jobs.py --pattern "case_pp*" --yes
```

## Common Workflows

### 1. Submit Everything

```bash
# Pure solid
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --validate --yes

# Porous flow
cd ../../pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --validate --yes
```

### 2. Submit Only Pure Solid Cases

```bash
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --validate --yes
```

### 3. Submit Only Porous Flow Cases

```bash
cd pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --validate --yes
```

### 4. Submit Specific Case Types

```bash
# Only pore pressure cases
cd pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --pattern "case_pp*" --validate --yes

# Only cf1 cases across both types
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --pattern "case_cf1*" --yes
cd ../../pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --pattern "case_cf1*" --yes
```

## Validation

Both scripts support comprehensive file validation:

```bash
# Validate pure solid cases
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --validate

# Validate porous flow cases
cd pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --validate
```

**What is validated:**
- Submit scripts exist (`.sh` files)
- Input files exist (`.i` files)
- Mesh files referenced in input files exist
- Files are in the correct locations

## Script Options

### Python Script (`submit_all_jobs.py`)

```bash
--validate, --check       # Validate all required files exist
--dry-run                 # Preview without submitting
--yes, -y                 # Skip confirmation prompt
--case <name>             # Submit specific case only
--pattern <pattern>       # Submit cases matching pattern
--delay <seconds>         # Delay between submissions
--script-name <name>      # Submit specific script name only
```

### Bash Script (`submit_all_jobs.sh`)

```bash
--validate, --check       # Validate all required files exist
--dry-run                 # Preview without submitting
--yes, -y                 # Skip confirmation prompt
--help, -h                # Show help message
<pattern>                 # Submit cases matching pattern (positional)
```

## Examples

### Example 1: Validate and submit all pure solid cases

```bash
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --validate && ./submit_all_jobs.py --yes
```

### Example 2: Validate and submit all porous flow cases

```bash
cd pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --validate && ./submit_all_jobs.py --yes
```

### Example 3: Submit only domain1x cases from both types

```bash
# Pure solid domain1x cases
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --pattern "*_domain1x" --yes

# Porous flow domain1x cases
cd ../../pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --pattern "*_domain1x" --yes
```

### Example 4: Quick validation check for both

```bash
# Check pure solid
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.sh --validate

# Check porous flow
cd ../../pf_code2d_porousflow/parametric_study
./submit_all_jobs.sh --validate
```

## Monitoring Jobs

After submitting jobs from either directory:

```bash
# Check all your jobs
squeue -u $USER

# Watch job queue (updates every 2 seconds)
watch -n 2 'squeue -u $USER'

# Cancel all jobs
scancel -u $USER

# Check output files (in respective directory)
tail -f case_*/*.o*
```

## Troubleshooting

### Wrong Directory

If you run the script and don't see the expected cases, check your current directory:

```bash
pwd  # Should show either pf_code2d_puresolid or pf_code2d_porousflow
```

### Missing Files

If validation fails, the script will tell you exactly which files are missing and in which case directory.

### No Cases Found

If the script says "No case folders found", you're likely in the wrong directory or the case directories haven't been created yet.

## Summary

- **Two independent sets of scripts**: One for pure solid, one for porous flow
- **Automatic directory detection**: Scripts only see files in their own directory
- **Same features in both**: Validation, dry-run, pattern matching, etc.
- **Navigate to the right directory** before running scripts
- **Use `--validate`** to check files before submitting
- **Use `--pattern`** to submit subsets of cases

For detailed documentation on all features, see `README_SUBMIT_JOBS.md` in each directory.
