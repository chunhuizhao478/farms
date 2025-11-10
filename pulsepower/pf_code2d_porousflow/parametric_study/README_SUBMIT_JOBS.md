# Job Submission Scripts

This directory contains scripts to submit all SLURM jobs for parametric studies.

## Available Scripts

### 1. Python Script: `submit_all_jobs.py`

Feature-rich Python script with comprehensive validation and control options.

**Features:**
- File validation to check if all required files exist
- Dry-run mode to preview submissions
- Pattern matching for selective submission
- Confirmation prompts (can be skipped)
- Delay between submissions
- Detailed error reporting

**Usage:**

```bash
# Validate all files before submission
./submit_all_jobs.py --validate

# Preview what would be submitted (dry run)
./submit_all_jobs.py --dry-run

# Submit all jobs (with confirmation)
./submit_all_jobs.py

# Validate and submit with confirmation
./submit_all_jobs.py --validate --yes

# Submit specific pattern
./submit_all_jobs.py --pattern "case_cf*"

# Submit specific case
./submit_all_jobs.py --case case_cf1_domain1x

# Submit without confirmation + 2 second delay
./submit_all_jobs.py --yes --delay 2
```

### 2. Bash Script: `submit_all_jobs.sh`

Simpler, faster bash alternative with validation support.

**Features:**
- File validation
- Dry-run mode
- Pattern matching
- Colorized output

**Usage:**

```bash
# Validate all files
./submit_all_jobs.sh --validate

# Dry run
./submit_all_jobs.sh --dry-run

# Submit all jobs
./submit_all_jobs.sh

# Submit with confirmation skip
./submit_all_jobs.sh --yes

# Submit specific pattern
./submit_all_jobs.sh "case_cf*"

# Submit specific case
./submit_all_jobs.sh case_cf1_domain1x
```

## File Validation

Both scripts can validate that all required files exist before submission:

**What is checked:**
- Submit scripts (`.sh` files)
- Input files (`.i` files)
- Mesh files (`.msh`, `.e`, `.exo`, `.mesh`)
- Referenced files in input files

**Example validation output:**

```
================================================================================
Validating required files...
================================================================================

  ✓ [case_cf1_domain1x] submit_elasticity.sh
  ✓ [case_cf1_domain2x] submit_elasticity.sh
  ✗ [case_cf1_domain5x] submit_elasticity.sh
      ✗ Missing Mesh file: fieldscale_test1_2d.msh
  ✓ [case_cf5_domain1x] submit_elasticity.sh

================================================================================
Validation Summary
================================================================================

Total jobs: 4
Valid: 3
Invalid (missing files): 1

✗ Validation FAILED - Some required files are missing
```

## Workflow

### Recommended Workflow

1. **Generate SLURM job scripts** (if not already done):
   ```bash
   ./generate_slurm_jobs.py
   ```

2. **Validate all files exist**:
   ```bash
   ./submit_all_jobs.py --validate
   ```

3. **Preview submissions** (optional):
   ```bash
   ./submit_all_jobs.py --dry-run
   ```

4. **Submit jobs**:
   ```bash
   ./submit_all_jobs.py --yes
   ```

### Quick Workflow

For experienced users who trust their setup:

```bash
# Validate and submit in one go
./submit_all_jobs.py --validate --yes
```

## Monitoring Jobs

After submission, monitor your jobs:

```bash
# Check job status
squeue -u $USER

# Watch job queue (updates every 2 seconds)
watch -n 2 'squeue -u $USER'

# Check detailed job info
scontrol show job <JOB_ID>

# Check output files
tail -f case_*/*.o*

# Cancel all jobs
scancel -u $USER

# Cancel specific job
scancel <JOB_ID>
```

## Common Issues

### Missing Files

If validation fails with missing files:

1. Check that mesh files are in the correct location
2. Verify file paths in input files are correct
3. Ensure files were generated if they should be auto-generated

### Permission Denied

If you get permission errors:

```bash
chmod +x submit_all_jobs.py
chmod +x submit_all_jobs.sh
```

### sbatch Command Not Found

This error means you're not on the HPC cluster. These scripts must be run on the cluster where SLURM is available.

## For Porous Flow Cases

To use these scripts with porous flow parametric studies:

```bash
cp submit_all_jobs.* ../pf_code2d_porousflow/parametric_study/
cd ../pf_code2d_porousflow/parametric_study/
./submit_all_jobs.py --validate
```

## Options Reference

### Python Script Options

| Option | Description |
|--------|-------------|
| `--validate`, `--check` | Validate all required files exist |
| `--dry-run` | Preview without submitting |
| `--yes`, `-y` | Skip confirmation prompt |
| `--case <name>` | Submit specific case only |
| `--pattern <pattern>` | Submit cases matching pattern |
| `--delay <seconds>` | Delay between submissions |
| `--script-name <name>` | Submit specific script name only |

### Bash Script Options

| Option | Description |
|--------|-------------|
| `--validate`, `--check` | Validate all required files exist |
| `--dry-run` | Preview without submitting |
| `--yes`, `-y` | Skip confirmation prompt |
| `--help`, `-h` | Show help message |
| `<pattern>` | Submit cases matching pattern (positional) |

## Exit Codes

- `0`: Success
- `1`: Error (missing files, submission failure, etc.)

## Examples

### Example 1: Validate before submitting all jobs

```bash
./submit_all_jobs.py --validate && ./submit_all_jobs.py --yes
```

### Example 2: Submit only compression cases with validation

```bash
./submit_all_jobs.py --validate --pattern "case_cf*" --yes
```

### Example 3: Submit with 5 second delays between jobs

```bash
./submit_all_jobs.py --delay 5 --yes
```

### Example 4: Quick validation check

```bash
./submit_all_jobs.sh --validate
```
