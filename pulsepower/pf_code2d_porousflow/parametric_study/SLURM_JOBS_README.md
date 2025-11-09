# SLURM Job Generation for Porousflow Parametric Studies

## Overview

Porousflow simulations require a **two-step workflow**:
1. **Step 1:** Run `static_solve.i` to compute static equilibrium
2. **Step 2:** Run `elasticity.i` for dynamic simulation (using static results)

## Complete Workflow

### Step 1: Static Solve (Equilibrium)

```bash
cd /path/to/pf_code2d_porousflow/parametric_study

# Generate SLURM scripts for static_solve.i
python3 generate_slurm_jobs.py --input-file static_solve.i

# Transfer to HPC
scp -r . username@frontera.tacc.utexas.edu:/scratch1/.../parametric_study

# Submit jobs on HPC
cd /scratch1/.../parametric_study
for dir in case_*/; do
    cd "$dir"
    sbatch submit_static_solve.sh
    cd ..
done

# Monitor jobs
squeue -u $USER

# Wait for completion...
```

### Step 2: Update Elasticity Files

After static solve completes, extract energy values and update `elasticity.i`:

```bash
# Download results from HPC
scp -r username@frontera:/.../parametric_study .

# Run update script locally
python3 run_static_solve_and_update.py --all --skip-simulation

# Verify updates
head -n 3 case_cf1_domain1x/elasticity.i
# Should show updated energy values:
# fluid_elastic_energy_total_static = 8.081664e-05
# solid_elastic_energy_total_static = 5.408951e-03
# full_input_energy_static = 5.489768e-03
```

### Step 3: Dynamic Simulation

```bash
# Generate SLURM scripts for elasticity.i
python3 generate_slurm_jobs.py --input-file elasticity.i

# Transfer updated files to HPC
scp -r . username@frontera.tacc.utexas.edu:/scratch1/.../parametric_study

# Submit dynamic simulation jobs
cd /scratch1/.../parametric_study
for dir in case_*/; do
    cd "$dir"
    sbatch submit_elasticity.sh
    cd ..
done
```

## Quick Commands

### Generate All Job Scripts

```bash
# Static solve jobs
python3 generate_slurm_jobs.py --input-file static_solve.i

# Dynamic simulation jobs
python3 generate_slurm_jobs.py --input-file elasticity.i

# Generate for specific cases
python3 generate_slurm_jobs.py --pattern "case_cf*" --input-file static_solve.i
```

### Custom SLURM Parameters

```bash
# Smaller jobs for static solve
python3 generate_slurm_jobs.py \
    --input-file static_solve.i \
    --nodes 4 \
    --ntasks 100 \
    --time 12:00:00

# Larger jobs for dynamic simulation
python3 generate_slurm_jobs.py \
    --input-file elasticity.i \
    --nodes 8 \
    --ntasks 200 \
    --time 48:00:00
```

## Configuration

Edit `generate_slurm_jobs.py`:

```python
# IMPORTANT: Update this path!
HPC_PROJECT_ROOT = "/scratch1/10024/zhaochun/projects/farms_cdms_11022025"

# Relative path to this study
RELATIVE_STUDY_PATH = "pulsepower/pf_code2d_porousflow/parametric_study"

# SLURM defaults
SLURM_DEFAULTS = {
    'partition': 'normal',
    'nodes': 8,
    'ntasks': 200,
    'time': '48:00:00',
    'account': 'EAR20006',
    'mail_user': 'chunhui3@illinois.edu',
}
```

## Generated Files

Each case will have **two** SLURM scripts:

```
case_cf1_domain1x/
├── static_solve.i
├── elasticity.i
├── fracture.i
├── submit_static_solve.sh   # ← For Step 1
└── submit_elasticity.sh     # ← For Step 3
```

## Command-Line Options

| Option | Description |
|--------|-------------|
| `--input-file FILE` | **Required:** `static_solve.i` or `elasticity.i` |
| `--case CASE` | Generate for specific case |
| `--pattern PATTERN` | Generate for matching cases |
| `--nodes N` | Number of nodes |
| `--ntasks N` | Total MPI tasks |
| `--time HH:MM:SS` | Wall time |

## Typical Resource Requirements

### Static Solve (static_solve.i)
- **Nodes:** 2-4
- **Tasks:** 50-100
- **Time:** 6-12 hours
- **Memory:** Moderate

### Dynamic Simulation (elasticity.i)
- **Nodes:** 8-16
- **Tasks:** 200-400
- **Time:** 24-48 hours
- **Memory:** High

## Monitoring Jobs

```bash
# Check all jobs
squeue -u $USER

# Check specific job
scontrol show job JOBID

# Watch output in real-time
tail -f case_cf1_domain1x/*.o*

# Check if static solve completed successfully
grep "Finished Executing" case_cf1_domain1x/*.o*
```

## Data Transfer

```bash
# Upload to HPC
scp -r parametric_study username@frontera:/scratch1/.../

# Download results (after static solve)
scp -r username@frontera:/.../parametric_study/*.csv .

# Download all results (after completion)
scp -r username@frontera:/.../parametric_study .
```

## Checklist

### Before Running Static Solve
- [ ] Generated cases with `generate_parametric_cases.py`
- [ ] Updated `HPC_PROJECT_ROOT` in `generate_slurm_jobs.py`
- [ ] Generated SLURM scripts: `--input-file static_solve.i`
- [ ] Transferred files to HPC
- [ ] Verified input files exist in each case folder

### After Static Solve Completes
- [ ] Downloaded CSV results from HPC
- [ ] Ran `run_static_solve_and_update.py --all --skip-simulation`
- [ ] Verified `elasticity.i` files are updated (check first 3 lines)
- [ ] Generated elasticity SLURM scripts: `--input-file elasticity.i`
- [ ] Transferred updated files back to HPC

### Before Running Dynamic Simulation
- [ ] Confirmed `elasticity.i` has correct energy values
- [ ] Verified `static_solve_out.e` checkpoint exists
- [ ] Generated elasticity SLURM scripts
- [ ] Ready to submit elasticity jobs

## Troubleshooting

### Static solve completes but no CSV
**Problem:** `static_solve_csv.csv` not found

**Solution:**
- Check MOOSE output for errors
- Verify Outputs block in `static_solve.i`
- Check CSV filename matches script expectation

### Elasticity.i still has old energy values
**Problem:** First 3 lines not updated after static solve

**Solution:**
```bash
# Re-run update script
python3 run_static_solve_and_update.py --all --skip-simulation

# Manually check CSV exists
ls case_cf1_domain1x/static_solve_csv.csv

# Debug extraction
python3 run_static_solve_and_update.py --case case_cf1_domain1x --skip-simulation --dry-run
```

### Dynamic simulation fails at startup
**Problem:** Missing checkpoint file

**Solution:**
- Ensure static_solve.i completed successfully
- Check for `static_solve_out_cp/` directory
- Verify checkpoint files exist

### Jobs pending too long
**Solution:**
```bash
# Check queue status
showq

# Check allocation
sbalance -u $USER

# Try different partition
python3 generate_slurm_jobs.py --partition development --time 2:00:00
```

## Example Workflow Script

Create a file `run_all.sh` on HPC:

```bash
#!/bin/bash

# Step 1: Submit all static solve jobs
echo "Submitting static solve jobs..."
for dir in case_*/; do
    cd "$dir"
    sbatch submit_static_solve.sh
    cd ..
done

# Wait and check completion (run manually after jobs finish)
# Step 2: After downloading and updating locally, submit elasticity jobs
# for dir in case_*/; do
#     cd "$dir"
#     sbatch submit_elasticity.sh
#     cd ..
# done
```

## Best Practices

1. **Test one case first** before submitting all jobs
2. **Monitor disk space** - simulations produce large output files
3. **Use checkpoint files** - enable restart if jobs fail
4. **Clean up old results** before re-running
5. **Document parameters** - keep notes on what's running
6. **Check allocations** - ensure sufficient hours available

## Support

For issues:
1. Check MOOSE documentation: https://mooseframework.inl.gov
2. Review HPC documentation (TACC Frontera)
3. Check job output/error files
4. Contact HPC support for cluster issues
