# Clamping Feature for Nonlocal Damage/Breakage

## Problem Statement

After nonlocal averaging using weighted spatial averaging, the damage (α) and breakage (B) values may fall slightly outside the physical bounds of [0, 1] due to:

1. **Numerical interpolation artifacts**
2. **Weighted averaging at boundaries** where neighbors have different values
3. **Boundary effects** near domain edges
4. **Quadrature point interpolation** errors

This can cause:
- ❌ Unphysical material properties (negative damage, damage > 1)
- ❌ Numerical instabilities in constitutive equations
- ❌ Convergence issues

## Solution: Automatic Clamping

### Implementation

Modified `DiffusedDamageBreakageMaterialSubAppNonlocal.C` to automatically clamp nonlocal values:

```cpp
// Store unclamped values for diagnostics
_alpha_before_clamp_mat[_qp] = _alpha_nonlocal_mat[_qp];
_B_before_clamp_mat[_qp] = _B_nonlocal_mat[_qp];

// Clamp to [0, 1]
_alpha_nonlocal_mat[_qp] = std::max(0.0, std::min(1.0, _alpha_nonlocal_mat[_qp]));
_B_nonlocal_mat[_qp] = std::max(0.0, std::min(1.0, _B_nonlocal_mat[_qp]));
```

**Location**: `src/materials/cdbm_lagrangian/DiffusedDamageBreakageMaterialSubAppNonlocal.C:332-339`

---

## Diagnostic Output Variables

To track when and where clamping occurs, the following variables are available:

### Sub App Output Variables

| Variable | Description | Purpose |
|----------|-------------|---------|
| `alpha_before_clamp` | Nonlocal damage BEFORE clamping | See unclamped values |
| `B_before_clamp` | Nonlocal breakage BEFORE clamping | See unclamped values |
| `alpha_nonlocal_sub` | Nonlocal damage AFTER clamping | Final clamped values |
| `B_nonlocal_sub` | Nonlocal breakage AFTER clamping | Final clamped values |

### How to Check if Clamping Occurred

#### Method 1: Visual Comparison in Paraview

```
1. Load: dynamic_solve_sub_eps1em7_cd10_cycle_sub_app0.e
2. Plot: alpha_before_clamp
3. Plot: alpha_nonlocal_sub
4. Calculator: abs(alpha_before_clamp - alpha_nonlocal_sub)
   → Non-zero values indicate clamping occurred
```

#### Method 2: Calculator Filter for Out-of-Bounds Detection

**Damage too high (α > 1):**
```python
Calculator: if(alpha_before_clamp > 1.0, alpha_before_clamp - 1.0, 0.0)
```

**Damage too low (α < 0):**
```python
Calculator: if(alpha_before_clamp < 0.0, -alpha_before_clamp, 0.0)
```

**Breakage too high (B > 1):**
```python
Calculator: if(B_before_clamp > 1.0, B_before_clamp - 1.0, 0.0)
```

**Breakage too low (B < 0):**
```python
Calculator: if(B_before_clamp < 0.0, -B_before_clamp, 0.0)
```

#### Method 3: Python Analysis

```python
import numpy as np

# Load data
alpha_before = data['alpha_before_clamp']
alpha_after = data['alpha_nonlocal_sub']

# Find where clamping occurred
clamped = np.abs(alpha_before - alpha_after) > 1e-10

# Statistics
n_clamped = np.sum(clamped)
pct_clamped = 100 * n_clamped / len(alpha_before)

print(f"Clamped points: {n_clamped} ({pct_clamped:.2f}%)")

# Where clamping occurred
too_high = alpha_before > 1.0
too_low = alpha_before < 0.0

print(f"Values > 1.0: {np.sum(too_high)}")
print(f"Values < 0.0: {np.sum(too_low)}")

# Maximum deviation
max_deviation = np.max(np.abs(alpha_before - alpha_after))
print(f"Max deviation: {max_deviation:.6e}")
```

---

## Expected Behavior

### ✅ Normal Case (Small Deviations)

Typical deviations are **very small** (< 0.01):

```
alpha_before_clamp = [0.0, 0.5, 0.98, 1.003, -0.001]
                              ↓ clamp ↓
alpha_nonlocal_sub = [0.0, 0.5, 0.98, 1.000,  0.000]
                              ↑ changed ↑
```

**This is fine!** Small numerical errors are expected and automatically corrected.

### ⚠️ Warning Case (Moderate Deviations)

If deviations are **moderate** (0.01 - 0.1):

```
alpha_before_clamp = [0.0, 0.5, 0.85, 1.08, -0.05]
                              ↓ clamp ↓
alpha_nonlocal_sub = [0.0, 0.5, 0.85, 1.00,  0.00]
                              ↑ significant ↑
```

**Action**: Check nonlocal averaging parameters:
- Reduce `radius` if too large
- Check `length_scale` is appropriate
- Verify boundary conditions

### ❌ Problem Case (Large Deviations)

If deviations are **large** (> 0.1):

```
alpha_before_clamp = [0.0, 0.5, 0.7, 1.5, -0.3]
                           ↓ BIG clamp ↓
alpha_nonlocal_sub = [0.0, 0.5, 0.7, 1.0,  0.0]
                           ↑ TOO MUCH ↑
```

**Action**: Something is wrong! Check:
1. Local damage evolution - is it producing valid [0,1] values?
2. Nonlocal averaging weights - are they normalized?
3. Material properties - are there NaN or Inf values?

---

## Visualization Workflow

### Step 1: Check if Clamping Occurred

```
Paraview → Calculator:
  Expression: abs(alpha_before_clamp - alpha_nonlocal_sub) > 1e-6
  Name: "was_clamped"

Color by: "was_clamped"
→ Blue (0) = no clamping
→ Red (1) = clamping occurred
```

### Step 2: Quantify the Clamping Effect

```
Paraview → Calculator:
  Expression: alpha_before_clamp - alpha_nonlocal_sub
  Name: "clamping_correction"

Color scale:
→ Green (≈0) = no correction needed
→ Red (>0) = clamped from above (was > 1.0)
→ Blue (<0) = clamped from below (was < 0.0)
```

### Step 3: Identify Problem Regions

```
Paraview → Threshold:
  Scalar: "clamping_correction"
  Range: [0.01, 1.0]  # Show only significant corrections

→ Highlights regions where clamping was significant
```

---

## When to Worry About Clamping

### ✅ Don't Worry If:

- Clamping affects < 1% of domain
- Deviations are < 0.01
- Only occurs at domain boundaries
- Temporary during transients

### ⚠️ Investigate If:

- Clamping affects > 10% of domain
- Deviations are 0.01 - 0.1
- Persistent throughout simulation
- Affects critical regions (damage zones)

### ❌ Stop and Debug If:

- Clamping affects > 50% of domain
- Deviations are > 0.1
- Occurs in low-damage regions
- Results look unphysical

---

## Output File Structure

```
dynamic_solve_sub_eps1em7_cd10_cycle_sub_app0.e
  ├─ alpha_damagedvar_sub      (local damage, evolved)
  ├─ alpha_local_output        (copy for comparison)
  ├─ alpha_before_clamp        ← BEFORE clamping (may be out of bounds)
  ├─ alpha_nonlocal_sub        ← AFTER clamping (guaranteed [0,1])
  ├─ B_damagedvar_sub          (local breakage, evolved)
  ├─ B_local_output            (copy for comparison)
  ├─ B_before_clamp            ← BEFORE clamping (may be out of bounds)
  └─ B_nonlocal_sub            ← AFTER clamping (guaranteed [0,1])
```

---

## Complete Comparison Pipeline

### Full Variable Chain:

```
DAMAGE:
  alpha_damagedvar_sub       ← [1] Local damage (sharp)
         ↓
  alpha_local_output         ← [2] Copy for easy viewing
         ↓
  [Nonlocal Averaging]
         ↓
  alpha_before_clamp         ← [3] Averaged (may be out of [0,1])
         ↓
  [Clamping to [0,1]]
         ↓
  alpha_nonlocal_sub         ← [4] Final smoothed & bounded
         ↓
  [Transfer to Main App]
         ↓
  alpha_damagedvar_aux       ← [5] Used in main app
```

```
BREAKAGE:
  B_damagedvar_sub           ← [1] Local breakage (sharp)
         ↓
  B_local_output             ← [2] Copy for easy viewing
         ↓
  [Nonlocal Averaging]
         ↓
  B_before_clamp             ← [3] Averaged (may be out of [0,1])
         ↓
  [Clamping to [0,1]]
         ↓
  B_nonlocal_sub             ← [4] Final smoothed & bounded
         ↓
  [Transfer to Main App]
         ↓
  B_damagedvar_aux           ← [5] Used in main app
```

---

## Quick Diagnostic Checklist

- [ ] Load sub app output file
- [ ] Plot `alpha_before_clamp` and `alpha_nonlocal_sub`
- [ ] Create Calculator: `abs(alpha_before_clamp - alpha_nonlocal_sub)`
- [ ] Check maximum difference (should be < 0.01)
- [ ] Check percentage of domain affected (should be < 1%)
- [ ] Verify clamping only at boundaries or extreme damage zones
- [ ] Same checks for breakage (B)

---

## Technical Notes

### Clamping Algorithm

```cpp
// Simple, robust clamping
clamped_value = std::max(0.0, std::min(1.0, unclamped_value));

// Equivalent to:
if (unclamped_value < 0.0)
  clamped_value = 0.0;
else if (unclamped_value > 1.0)
  clamped_value = 1.0;
else
  clamped_value = unclamped_value;
```

### Why Clamping is Needed

Nonlocal averaging uses weighted sums:

```
α_nonlocal(x) = Σ w(x,x') * α(x') / Σ w(x,x')
```

Even if all `α(x') ∈ [0,1]`, numerical issues can cause:
- Interpolation at quadrature points → small errors
- Floating point arithmetic → roundoff errors
- Boundary weights not exactly normalized → extrapolation

Clamping ensures physical validity without affecting results significantly.

---

## References

- Material file: `DiffusedDamageBreakageMaterialSubAppNonlocal.C:326-339`
- Input file: `dynamic_solve_sub_eps1em7_cd10_cycle.i:149-157, 276-288`
- Comparison guide: `COMPARISON_VARIABLES_GUIDE.md`

---

*Feature added to ensure numerical stability and physical validity of nonlocal averaging results.*
