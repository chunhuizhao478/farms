# Quick Reference: Output Variables for Comparison

## Sub App Output File: `dynamic_solve_sub_eps1em7_cd10_cycle_sub_app0.e`

### Damage Variables (α):

```
┌─────────────────────────────────────────────────────────────┐
│                    DAMAGE COMPARISON                         │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  BEFORE AVERAGING (Local):                                  │
│    ✓ alpha_damagedvar_sub      ← Primary evolution variable │
│    ✓ alpha_local_output        ← Copy for easy viewing      │
│                                                              │
│                         ↓↓↓                                  │
│                  [Nonlocal Averaging]                        │
│                  • length_scale = 200                        │
│                  • radius = 400                              │
│                  • weights = BAZANT                          │
│                         ↓↓↓                                  │
│                                                              │
│  AFTER AVERAGING (Before Clamping):                         │
│    ✓ alpha_before_clamp        ← May be out of [0,1] bounds │
│                                                              │
│                         ↓↓↓                                  │
│                     [Clamping to [0,1]]                      │
│                         ↓↓↓                                  │
│                                                              │
│  AFTER CLAMPING (Nonlocal):                                 │
│    ✓ alpha_nonlocal_sub        ← Final smoothed & bounded   │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Breakage Variables (B):

```
┌─────────────────────────────────────────────────────────────┐
│                   BREAKAGE COMPARISON                        │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  BEFORE AVERAGING (Local):                                  │
│    ✓ B_damagedvar_sub          ← Primary evolution variable │
│    ✓ B_local_output            ← Copy for easy viewing      │
│                                                              │
│                         ↓↓↓                                  │
│                  [Nonlocal Averaging]                        │
│                  • length_scale = 200                        │
│                  • radius = 400                              │
│                  • weights = BAZANT                          │
│                         ↓↓↓                                  │
│                                                              │
│  AFTER AVERAGING (Before Clamping):                         │
│    ✓ B_before_clamp            ← May be out of [0,1] bounds │
│                                                              │
│                         ↓↓↓                                  │
│                     [Clamping to [0,1]]                      │
│                         ↓↓↓                                  │
│                                                              │
│  AFTER CLAMPING (Nonlocal):                                 │
│    ✓ B_nonlocal_sub            ← Final smoothed & bounded   │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## Main App Output File: `unified_solve_main_eps1em7_cd10_cycle_out.e`

### Received Nonlocal Variables:

```
┌─────────────────────────────────────────────────────────────┐
│              MAIN APP (Uses Nonlocal Values)                 │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ✓ alpha_damagedvar_aux  ← NONLOCAL damage from sub app    │
│                             (should match alpha_nonlocal_sub)│
│                                                              │
│  ✓ B_damagedvar_aux      ← NONLOCAL breakage from sub app  │
│                             (should match B_nonlocal_sub)   │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## Visualization Workflow in Paraview

### Step 1: Load Sub App Output
```
File → Open → dynamic_solve_sub_eps1em7_cd10_cycle_sub_app0.e
```

### Step 2: Create Comparison Plots

#### Option A: Side-by-Side Comparison
1. Split view horizontally
2. **Left panel**: Show `alpha_local_output` (local)
3. **Right panel**: Show `alpha_nonlocal_sub` (nonlocal)
4. Use same color scale for both

#### Option B: Difference Plot
1. Calculator filter: `alpha_nonlocal_sub - alpha_local_output`
2. This shows the smoothing effect
3. Positive values = averaging increased damage
4. Negative values = averaging decreased damage

### Step 3: Verify Transfer to Main App
1. Load main app output: `unified_solve_main_eps1em7_cd10_cycle_out.e`
2. Plot `alpha_damagedvar_aux`
3. Compare visually with `alpha_nonlocal_sub` from sub app
4. They should be identical (within numerical precision)

---

## Expected Observations

### ✅ What You Should See:

1. **Local fields are sharper**:
   - `alpha_local_output` and `B_local_output` have steep gradients
   - Localization bands are narrow
   - High spatial frequency content

2. **Nonlocal fields are smoother**:
   - `alpha_nonlocal_sub` and `B_nonlocal_sub` have gentle gradients
   - Localization bands are wider (influenced by `length_scale = 200`)
   - Spatial smoothing within `radius = 400`

3. **Peak values may differ**:
   - Local peaks can be higher (concentrated)
   - Nonlocal peaks are lower but wider (distributed)
   - Integrated quantities should be similar

4. **Mesh sensitivity reduced**:
   - Nonlocal results should be less sensitive to mesh refinement
   - Local results may show mesh-dependent localization

### ❌ What Would Indicate Problems:

1. `alpha_nonlocal_sub` identical to `alpha_local_output` → averaging not working
2. `alpha_damagedvar_aux` ≠ `alpha_nonlocal_sub` → transfer failed
3. Excessive smoothing → `radius` too large
4. No visible smoothing → `radius` too small or averaging not applied

---

## Quantitative Comparison (Python/Paraview)

### Calculate Smoothing Metrics:

```python
import numpy as np

# Load sub app data
alpha_local = data['alpha_local_output']
alpha_nonlocal = data['alpha_nonlocal_sub']

# 1. Maximum difference
max_diff = np.max(np.abs(alpha_nonlocal - alpha_local))
print(f"Max smoothing effect: {max_diff}")

# 2. RMS difference
rms_diff = np.sqrt(np.mean((alpha_nonlocal - alpha_local)**2))
print(f"RMS smoothing effect: {rms_diff}")

# 3. Spatial gradient reduction
grad_local = np.gradient(alpha_local)
grad_nonlocal = np.gradient(alpha_nonlocal)
grad_ratio = np.linalg.norm(grad_nonlocal) / np.linalg.norm(grad_local)
print(f"Gradient reduction: {grad_ratio:.3f}")

# 4. Conservation check (should be close to 1.0)
integral_ratio = np.sum(alpha_nonlocal) / np.sum(alpha_local)
print(f"Mass conservation: {integral_ratio:.6f}")
```

---

## Output File Locations

```
/Users/chunhuizhao/projects/farms_cdbm_implicit/
  └── development/paper_draft_slip_spectrum/cases/cycletest_dev/
      ├── unified_solve_main_eps1em7_cd10_cycle_out.e        ← Main app
      ├── unified_solve_main_eps1em7_cd10_cycle_out_sub_app0.e  ← Sub app
      └── COMPARISON_VARIABLES_GUIDE.md  ← This file
```

---

## Variable Summary Table

| Variable | Type | Location | Description | When to Use |
|----------|------|----------|-------------|-------------|
| `alpha_damagedvar_sub` | Variable | Sub app | Local damage (evolved) | Primary evolution |
| `alpha_local_output` | AuxVar | Sub app | Copy of local damage | **Easy comparison** |
| `alpha_before_clamp` | AuxVar | Sub app | Averaged BEFORE clamping | **Check for out-of-bounds** |
| `alpha_nonlocal_sub` | AuxVar | Sub app | Averaged AFTER clamping | **Final smoothed result** |
| `B_damagedvar_sub` | Variable | Sub app | Local breakage (evolved) | Primary evolution |
| `B_local_output` | AuxVar | Sub app | Copy of local breakage | **Easy comparison** |
| `B_before_clamp` | AuxVar | Sub app | Averaged BEFORE clamping | **Check for out-of-bounds** |
| `B_nonlocal_sub` | AuxVar | Sub app | Averaged AFTER clamping | **Final smoothed result** |
| `alpha_damagedvar_aux` | AuxVar | Main app | Nonlocal damage (received) | Main app constitutive |
| `B_damagedvar_aux` | AuxVar | Main app | Nonlocal breakage (received) | Main app constitutive |

---

## Quick Checks

### ✓ Is nonlocal averaging working?
→ Plot `alpha_local_output` vs `alpha_nonlocal_sub` side-by-side
→ Nonlocal should be smoother

### ✓ Is transfer working?
→ Compare `alpha_nonlocal_sub` (sub app) vs `alpha_damagedvar_aux` (main app)
→ Should be identical

### ✓ Are results reasonable?
→ Check peak values: nonlocal peaks should be lower but wider
→ Check integrated values: should be similar
→ Check spatial extent: nonlocal should affect ≈ 2×radius around damage

---

*Quick tip: Start by visualizing `alpha_local_output` vs `alpha_nonlocal_sub` to immediately see the nonlocal smoothing effect!*
