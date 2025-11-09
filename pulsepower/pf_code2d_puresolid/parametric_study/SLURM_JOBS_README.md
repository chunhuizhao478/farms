# SLURM Job Generation for Pure Solid Parametric Studies

## Quick Start

### 1. Generate SLURM Job Scripts

```bash
cd /path/to/pf_code2d_puresolid/parametric_study

# Generate for all cases (default: elasticity.i)
python3 generate_slurm_jobs.py

# Generate for specific cases
python3 generate_slurm_jobs.py --pattern "case_cf*"
python3 generate_slurm_jobs.py --case case_cf1_domain1x

# Custom SLURM parameters
python3 generate_slurm_jobs.py --nodes 4 --ntasks 100 --time 24:00:00
```

### 2. Transfer to HPC

```bash
# Update HPC_PROJECT_ROOT in generate_slurm_jobs.py first!
# Then transfer to cluster
scp -r . username@frontera.tacc.utexas.edu:/scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_puresolid/parametric_study
```

### 3. Submit Jobs on HPC

```bash
# SSH to cluster
ssh username@frontera.tacc.utexas.edu

# Navigate to study directory
cd /scratch1/10024/zhaochun/projects/farms_cdms_11022025/pulsepower/pf_code2d_puresolid/parametric_study

# Submit single job
cd case_cf1_domain1x
sbatch submit_elasticity.sh

# Submit all jobs
cd ..
for dir in case_*/; do
    cd "$dir"
    sbatch submit_elasticity.sh
    cd ..
done
```

### 4. Monitor Jobs

```bash
# Check job status
squeue -u $USER

# Check detailed job info
scontrol show job JOBID

# Check output
tail -f case_cf1_domain1x/*.o*
```

## Configuration

Edit the **CONFIGURATION** section in `generate_slurm_jobs.py`:

```python
# HPC project root path (IMPORTANT: Update this!)
HPC_PROJECT_ROOT = "/scratch1/10024/zhaochun/projects/farms_cdms_11022025"

# SLURM defaults
SLURM_DEFAULTS = {
    'partition': 'normal',      # Queue name
    'nodes': 8,                  # Number of nodes
    'ntasks': 200,               # Total MPI tasks
    'time': '48:00:00',         # Wall time
    'account': 'EAR20006',      # Allocation
    'mail_user': 'chunhui3@illinois.edu',
}
```

## Command-Line Options

| Option | Description |
|--------|-------------|
| `--case CASE` | Generate for specific case |
| `--pattern PATTERN` | Generate for cases matching pattern |
| `--nodes N` | Number of nodes (default: 8) |
| `--ntasks N` | Total MPI tasks (default: 200) |
| `--time HH:MM:SS` | Wall time (default: 48:00:00) |
| `--partition NAME` | Queue name (default: normal) |
| `--account NAME` | Project allocation |
| `--input-file FILE` | Input file to run (default: elasticity.i) |

## Generated Files

Each case folder will contain:
```
case_cf1_domain1x/
├── elasticity.i
├── fracture.i
└── submit_elasticity.sh  # ← Generated SLURM script
```

## Example SLURM Script

```bash
#!/bin/bash
#SBATCH -J case_cf1_domain1x_elasticity
#SBATCH -o case_cf1_domain1x_elasticity.o%j
#SBATCH -e case_cf1_domain1x_elasticity.e%j
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 200
#SBATCH -t 48:00:00
#SBATCH --mail-type=all
#SBATCH -A EAR20006
#SBATCH --mail-user=chunhui3@illinois.edu

# Load modules
module swap intel gcc
module load cuda
export CC=mpicc CXX=mpicxx FC=mpif90

# Run simulation
ibrun ./farms-opt -i /path/to/case_cf1_domain1x/elasticity.i --allow-unused
```

## Tips

1. **Update HPC_PROJECT_ROOT** before generating scripts
2. **Test one case** before submitting all jobs
3. **Monitor queue limits** - don't submit too many jobs at once
4. **Check module availability** on your HPC system
5. **Verify paths** in generated scripts before submission

## Troubleshooting

### Job fails immediately
- Check module commands are correct for your HPC
- Verify executable path
- Check input file exists

### Jobs stuck in queue
- Check queue limits: `squeue -u $USER`
- Verify account has allocation: `sbalance -u $USER`
- Check partition availability

### Out of memory errors
- Increase nodes: `--nodes 16`
- Decrease tasks per node

### Wrong paths in script
- Update `HPC_PROJECT_ROOT` in `generate_slurm_jobs.py`
- Regenerate scripts: `python3 generate_slurm_jobs.py`
