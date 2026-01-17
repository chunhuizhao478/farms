# Large Deformation Convergence Fixes

## Problem
Negative Jacobian errors causing element inversion during large deformation analysis.
- Original failure at t=7.9, element 3005643
- After first fixes: failure at t=9.13, element 1796752
- Location: near loading point (x≈0.014, y≈0.0007, z≈0.0002)

## Applied Fixes (Most Aggressive)

### 1. **Reduced Loading Rate by 10x** ✓
```
Before: expression = '-1e-4 * t'
After:  expression = '-1e-5 * t'
```
**Reason**: Loading was too fast for large deformation. Elements need time to redistribute stress.

### 2. **More Conservative Time Stepping** ✓
```
dt (initial):        0.01  → 0.001 (10x smaller)
optimal_iterations:  8    → 6 (more conservative)
growth_factor:       1.5  → 1.2 (slower growth)
dtmin:              1e-6  → 1e-7 (can go smaller)
dtmax:              0.1   → 0.05 (limit max step)
```

### 3. **Switched Decomposition Method** ✓
```
Before: decomposition = spectral
After:  decomposition = VOLDEV
```
**Reason**: Spectral decomposition requires eigenvalue decomposition which can fail for distorted elements. VOLDEV is more robust.

### 4. **Extended Simulation Time** ✓
```
Before: end_time = 30
After:  end_time = 300
```
**Reason**: Reduced loading rate means we need 10x more time to reach same displacement.

### 5. **Explicit Initial Conditions** ✓
```
Added: initial_condition = 0.0 to all displacement variables
```
**Reason**: Ensure clean start with zero displacement.

### 6. **Added Checkpoint/Restart** ✓
```
checkpoint = true
num_checkpoint_files = 2
```
**Reason**: Can restart from last good time step if failure occurs.

### 7. **Using Coarse Mesh** ✓
```
file = '../../meshfile/mesh_wohole_3d_quicktest.msh'
```
**Reason**: Faster iteration for testing convergence parameters.

---

## Summary of Current Settings

| Parameter | Value | Purpose |
|-----------|-------|---------|
| Loading rate | 1e-5 m/s | 10x slower |
| Initial dt | 0.001 | Very small start |
| Min dt | 1e-7 | Can cut back far |
| Max dt | 0.05 | Limited growth |
| Decomposition | VOLDEV | Robust method |
| eta (residual) | 1e-3 | Prevent full degradation |
| Line search | basic | Prevent overshoot |
| Mesh | coarse | Fast testing |

---

## What to Expect

### If It Converges:
- Should reach t ~ 30-50 before damage initiates
- At t=30: displacement = -0.0003 m (much less than before)
- Adaptive dt will gradually increase as solution stabilizes
- Check max(d) to see when fracture starts

### If It Still Fails:
Try these additional measures in order:

#### Option A: Even Slower Loading
```
expression = '-1e-6 * t'  # 100x slower than original
end_time = 3000
```

#### Option B: Switch to No Decomposition
```
decomposition = NONE
```
Even simpler, no energy split at all.

#### Option C: Test Without Phase Field
Comment out MultiApp temporarily:
```
# [MultiApps]
#   [fracture]
#   ...
# []
```
This tests if the issue is purely mechanical.

#### Option D: Use Direct Solver
```
petsc_options_iname = '-pc_type -snes_type'
petsc_options_value = 'lu vinewtonrsls'
```
More robust but slower.

#### Option E: Increase Residual Stiffness
```
parameter_values = '2 0.01'  # eta = 0.01 (1% residual)
```

---

## Running the Simulation

### Clean Start
```bash
# Remove old output files
rm -f elasticity_out*.e elasticity_out*.cpr

# Run simulation
mpiexec -n 4 <your-executable> -i elasticity.i
```

### Monitor Progress
```bash
# Watch the output
tail -f <output_log>

# Look for:
# - Time step size (should start at 0.001)
# - Nonlinear iterations (should be 5-8)
# - Any cutbacks (indicates difficulty)
```

### Restart from Checkpoint (if needed)
```bash
# If simulation fails at t=50, restart from last checkpoint
mpiexec -n 4 <your-executable> -i elasticity.i --recover
```

---

## Diagnostics After Run

### Check Maximum Displacement
```bash
# In ParaView or using Python
# Look for max(disp_y)
# Should be much smaller than original run
```

### Check Damage Initiation
```bash
# Look for max(d) > 0
# Note the time when fracture starts
```

### Check Time Step Evolution
```bash
# From output file
grep "Time Step" <log_file>
# Should see dt gradually increase
# Cutbacks indicate convergence difficulty
```

### Analyze Failed Element Location
If still failing:
1. Load mesh in ParaView
2. Search for failed element ID
3. Check mesh quality in that region
4. Consider local refinement

---

## Transition Back to Fine Mesh

Once coarse mesh runs successfully:

1. **Verify full simulation completes**
   ```bash
   # Let it run to end_time
   ```

2. **Check results make physical sense**
   - Displacement field smooth
   - Damage initiates at expected location
   - Force-displacement curve reasonable

3. **Switch to fine mesh**
   ```
   file = '../../meshfile/mesh_wohole_3d.msh'
   ```

4. **Start very conservatively**
   ```
   dt = 0.0001  # Even smaller for fine mesh
   dtmax = 0.01
   ```

5. **Consider using coarse mesh solution as initial condition**
   ```
   # Advanced: interpolate coarse solution to fine mesh
   ```

---

## Theory: Why These Fixes Help

### Negative Jacobian = Element Inversion
When deformation is too large, element can fold inside-out:
```
Original:     Inverted:
  3----4        4----3
  |    |        |    |
  1----2        2----1
```

### Prevention Strategies:

1. **Smaller time steps** → Smaller incremental deformation per step
2. **Slower loading** → Less total deformation at any given time
3. **Simpler decomposition** → Less complex tensor operations
4. **Residual stiffness** → Damaged elements maintain some strength
5. **Line search** → Prevents Newton from taking too large steps
6. **Adaptive dt** → Automatically reduces step when convergence is hard

### Energy Balance
In large deformation:
```
ψ = ψ_elastic + ψ_fracture
```
Must ensure:
- Elastic energy computed correctly (Hencky strain)
- Degradation doesn't cause complete loss of stiffness
- Energy split (VOLDEV) is robust to distortion

---

## Alternative Approach: Quasi-Static Ramp

If problems persist, consider ramped loading:

```python
[Functions]
  [func_loading]
    type = PiecewiseLinear
    x = '0   10   end_time'
    y = '0  -0.001  -0.003'  # Ramp up gradually
  []
[]
```

This gives elements time to adjust at the beginning.

---

## Contact Information

If issues persist after all these fixes:
1. Check mesh quality in the problematic region
2. Consider if boundary conditions might be causing local stress concentration
3. Verify material parameters are reasonable for large deformation
4. Consider if the problem physics requires plasticity (not just phase field)

---

## Quick Reference: Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| Fails immediately | Loading too fast | Reduce dt, reduce loading rate |
| Fails at same time | Mesh/BC issue | Check problematic element location |
| Fails randomly | Numerical precision | Use direct solver, tighter tolerances |
| Converges but wrong results | Wrong formulation | Check material model, BCs |
| Never converges | Problem too stiff | Increase eta, use NONE decomposition |

---

## Current Status

✓ Applied all aggressive fixes
✓ Using coarse mesh for testing
✓ Checkpoint enabled for restart
⏳ Ready to test

**Next Step**: Run and monitor convergence behavior.
