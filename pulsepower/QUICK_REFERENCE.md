# Quick Reference: Job Submission Scripts

## Files Created

### Pure Solid Directory
`pf_code2d_puresolid/parametric_study/`
- `submit_all_jobs.py` - Python submission script (9 cases)
- `submit_all_jobs.sh` - Bash submission script (9 cases)
- `README_SUBMIT_JOBS.md` - Detailed documentation

### Porous Flow Directory
`pf_code2d_porousflow/parametric_study/`
- `submit_all_jobs.py` - Python submission script (18 cases, 27 jobs)
- `submit_all_jobs.sh` - Bash submission script (18 cases, 27 jobs)
- `README_SUBMIT_JOBS.md` - Detailed documentation

### Guide Documents
- `pulsepower/SUBMIT_JOBS_GUIDE.md` - Complete usage guide
- `pulsepower/QUICK_REFERENCE.md` - This file

## Key Features

✓ **Isolated Operation**: Each script only sees files in its own directory
✓ **File Validation**: Check all required files exist before submission
✓ **Dry Run Mode**: Preview what would be submitted
✓ **Pattern Matching**: Submit specific subsets of cases
✓ **Automatic Detection**: Scripts auto-detect their location

## Common Commands

### Pure Solid Cases (9 cases)

```bash
cd pf_code2d_puresolid/parametric_study

# Validate files
./submit_all_jobs.py --validate

# Preview submission
./submit_all_jobs.py --dry-run

# Submit all
./submit_all_jobs.py --yes

# Submit specific pattern
./submit_all_jobs.py --pattern "case_cf5*" --yes
```

### Porous Flow Cases (18 cases)

```bash
cd pf_code2d_porousflow/parametric_study

# Validate files
./submit_all_jobs.py --validate

# Preview submission
./submit_all_jobs.py --dry-run

# Submit all
./submit_all_jobs.py --yes

# Submit only compression factor cases
./submit_all_jobs.py --pattern "case_cf*" --yes

# Submit only pore pressure cases
./submit_all_jobs.py --pattern "case_pp*" --yes
```

## Verification

**Pure Solid**: Finds 9 jobs
```bash
cd pf_code2d_puresolid/parametric_study
./submit_all_jobs.py --dry-run
# Output: Found 9 job(s) to submit
```

**Porous Flow**: Finds 27 jobs (18 cases, some with multiple job scripts)
```bash
cd pf_code2d_porousflow/parametric_study
./submit_all_jobs.py --dry-run
# Output: Found 27 job(s) to submit
```

## Best Practice Workflow

1. **Navigate to desired directory**
   ```bash
   cd pf_code2d_[puresolid|porousflow]/parametric_study
   ```

2. **Validate all files exist**
   ```bash
   ./submit_all_jobs.py --validate
   ```

3. **Preview what will be submitted**
   ```bash
   ./submit_all_jobs.py --dry-run
   ```

4. **Submit jobs**
   ```bash
   ./submit_all_jobs.py --yes
   ```

5. **Monitor jobs**
   ```bash
   squeue -u $USER
   watch -n 2 'squeue -u $USER'
   ```

## Options Quick Reference

| Option | Python | Bash | Description |
|--------|--------|------|-------------|
| Validate files | `--validate` | `--validate` | Check all files exist |
| Dry run | `--dry-run` | `--dry-run` | Preview only |
| Skip confirm | `--yes` | `--yes` | No prompt |
| Pattern match | `--pattern "case_*"` | `case_*` | Filter cases |
| Specific case | `--case case_name` | N/A | One case only |
| Delay | `--delay 2` | N/A | Wait between jobs |
| Help | `--help` | `--help` | Show help |

## Troubleshooting

**No cases found?**
→ Check you're in the right directory with `pwd`

**Wrong number of cases?**
→ Make sure you're in the correct directory (puresolid vs porousflow)

**Validation fails?**
→ Script will show exactly which files are missing

**Want to test?**
→ Use `--dry-run` to preview without submitting

## Summary

- **2 independent directories** with identical script functionality
- **Pure solid**: 9 cases (cf1, cf5, cf10 across 3 domain sizes)
- **Porous flow**: 18 cases (cf + pp across 3 domain sizes)
- **Navigate first**, then run scripts
- **Always validate** before submitting to avoid job failures

For complete documentation, see:
- `README_SUBMIT_JOBS.md` in each directory
- `SUBMIT_JOBS_GUIDE.md` in pulsepower directory
