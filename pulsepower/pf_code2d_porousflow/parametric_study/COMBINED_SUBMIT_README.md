# Combined Static + Dynamic Submission Scripts

## Overview

New combined submission scripts have been created for all porous flow cases. Each `submit_combined.sh` script runs **both** simulation steps in a single SLURM job:

1. **Static solve** (`static_solve.i`) - runs first
2. **Dynamic solve** (`elasticity.i`) - runs after static solve completes

## Files Created

**Generator script:**
- `generate_combined_submit.py` - Generates combined submit scripts

**Combined submit scripts (18 total):**
- `case_*/submit_combined.sh` - One per case directory

## Advantages of Combined Scripts

✓ **Single job submission** - No need to wait for static solve to finish before submitting dynamic
✓ **Automatic sequencing** - Dynamic solve starts immediately after static completes
✓ **Error handling** - Job aborts if static solve fails
✓ **Single allocation** - Both steps use the same compute nodes
✓ **Simpler workflow** - Submit once, get both results

## Comparison: Separate vs Combined

### Separate Jobs (Old Way)

```bash
# Step 1: Submit static solve
sbatch submit_static_solve.sh

# Step 2: Wait for static to finish, then submit dynamic
# (manual intervention required)
sbatch submit_elasticity.sh
```

**Issues:**
- Need to monitor static solve completion
- Manual submission of second job
- Potential queue wait between jobs

### Combined Job (New Way)

```bash
# Submit once - runs both automatically
sbatch submit_combined.sh
```

**Benefits:**
- Single submission
- Automatic progression
- No queue wait between steps

## Usage on HPC Cluster

### Submit All Combined Jobs

```bash
cd /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_porousflow/parametric_study

# Using Python script (recommended)
./submit_all_jobs.py --script-name submit_combined.sh --yes

# Or using bash loop
for dir in case_*/; do
  if [ -f "$dir/submit_combined.sh" ]; then
    cd "$dir" && sbatch submit_combined.sh && cd ..
  fi
done
```

### Submit Specific Cases

```bash
# Submit only compression factor cases
./submit_all_jobs.py --pattern "case_cf*" --script-name submit_combined.sh --yes

# Submit only pore pressure cases
./submit_all_jobs.py --pattern "case_pp*" --script-name submit_combined.sh --yes

# Submit single case
cd case_cf1_domain1x
sbatch submit_combined.sh
```

## What Each Combined Script Does

Each `submit_combined.sh` performs these steps:

1. **Setup**
   - Load modules
   - Set environment variables
   - Change to case directory

2. **Static Solve**
   - Print start banner
   - Run: `ibrun ./farms-opt -i static_solve.i`
   - Check exit code
   - Verify output file exists (`static_solve_out.e`)
   - **Abort if failed**

3. **Dynamic Solve**
   - Print start banner
   - Run: `ibrun ./farms-opt -i elasticity.i`
   - Check exit code
   - **Abort if failed**

4. **Completion**
   - Print success banner
   - Exit with status code

## Output Files

Each combined job produces:

- `case_NAME_combined.oJOBID` - Standard output (both steps)
- `case_NAME_combined.eJOBID` - Standard error (both steps)
- `static_solve_out.e` - Static solve results
- `elasticity_out.e` - Dynamic solve results

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Watch all your jobs
watch -n 2 'squeue -u $USER'

# Check combined job output (live)
tail -f case_*/*_combined.o*

# Check specific case
tail -f case_cf1_domain1x/case_cf1_domain1x_combined.o*

# Check for errors
grep -i error case_*/*_combined.e*
```

## Regenerating Combined Scripts

If you need to update parameters or regenerate scripts:

```bash
# Edit configuration in generate_combined_submit.py
# Then run:
python3 generate_combined_submit.py
```

This will recreate all `submit_combined.sh` files with updated settings.

## Resource Usage

**Time allocation:** 48:00:00 (same as individual jobs)
- Static solve typically takes: ~2-4 hours
- Dynamic solve typically takes: ~40-44 hours
- Total fits within 48-hour limit

**Compute resources:**
- Nodes: 8
- Tasks: 200
- Partition: normal

## Validation

Before submitting to cluster, validate locally:

```bash
# Dry run to check all scripts exist
./submit_all_jobs.py --script-name submit_combined.sh --dry-run

# Validate files exist
./submit_all_jobs.py --script-name submit_combined.sh --validate
```

## Troubleshooting

### Static solve fails

Check the output file:
```bash
tail -100 case_NAME/case_NAME_combined.o*
```

Look for error messages after "STEP 1: Running static solve..."

### Dynamic solve fails

Check if static output exists:
```bash
ls -lh case_NAME/static_solve_out.e
```

Check output after "STEP 2: Running dynamic solve..."

### Job aborted unexpectedly

Check error file:
```bash
cat case_NAME/case_NAME_combined.e*
```

Look for SLURM errors or resource issues

## Comparison Table

| Feature | Separate Scripts | Combined Script |
|---------|-----------------|-----------------|
| Submissions needed | 2 per case | 1 per case |
| Manual intervention | Yes (between jobs) | No |
| Queue wait | 2× (once per job) | 1× (single job) |
| Error handling | Manual check | Automatic abort |
| Total jobs (18 cases) | 36 jobs | 18 jobs |
| Job monitoring | Track 2 job IDs | Track 1 job ID |
| Node allocation | Potentially different | Same nodes |

## Current Script Inventory

For each of the 18 porous flow cases, you now have:

- `submit_static_solve.sh` - Static only (old)
- `submit_elasticity.sh` - Dynamic only (old)
- `submit_combined.sh` - Both steps (new) **← Use this**

## Recommendation

**Use `submit_combined.sh` for all new runs** unless you specifically need to run only one step.

The old separate scripts are kept for backwards compatibility and special cases.
