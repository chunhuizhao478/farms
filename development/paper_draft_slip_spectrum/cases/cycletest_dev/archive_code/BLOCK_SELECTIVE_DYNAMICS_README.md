# Block-Selective Dynamic/Quasi-Dynamic Formulation

## Overview

This implementation uses **spatially selective inertia** to prevent boundary loading waves from polluting the interior domain results during dynamic rupture simulations.

## Problem Statement

### The Challenge

When applying boundary loading (displacement or stress) to simulate tectonic forcing:
- Boundary conditions generate **loading waves** that propagate into the domain
- These waves can interact with and contaminate the rupture dynamics
- In full dynamic mode, these boundary waves carry artificial inertial effects

### Solution

Implement **region-dependent formulation**:
- **Outer block (block 1)**: Always quasi-dynamic (no inertia)
- **Inner region (blocks 2, 3)**: Adaptive switching (dynamic ↔ quasi-dynamic)

This creates a **buffer zone** that absorbs boundary artifacts while allowing natural rupture dynamics in the fault region.

---

## Implementation

### Block Structure

```
┌─────────────────────────────────────────────────────┐
│  Block 1: Outer Region (30 km × 30 km)             │
│  - Always quasi-dynamic (radiation damping)         │
│  - No inertia effects                               │
│  - Sochacki sponge absorbing boundaries             │
│  ┌───────────────────────────────────────────────┐  │
│  │ Blocks 2, 3: Inner Region (20 km × 20 km)    │  │
│  │ - Adaptive dynamic/quasi-dynamic switching    │  │
│  │ - Inertia enabled when ε̇ > 1e-5              │  │
│  │ - Contains fault zone and damage region       │  │
│  │                                               │  │
│  │     [Damage Box: -15 km to +15 km]           │  │
│  │                                               │  │
│  └───────────────────────────────────────────────┘  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Kernel Configuration

#### 1. Inertia (Dynamic Mode) - Inner Region Only

```
[Kernels]
  [inertia_x]
    type = InertialForce
    variable = disp_x
    acceleration = accel_x
    velocity = vel_x
    beta = 0.25
    gamma = 0.5
    eta = 0
    enable = false
    block = '2 3'  # RESTRICTED TO INNER REGION
  []
  [inertia_y]
    type = InertialForce
    variable = disp_y
    acceleration = accel_y
    velocity = vel_y
    beta = 0.25
    gamma = 0.5
    eta = 0
    enable = false
    block = '2 3'  # RESTRICTED TO INNER REGION
  []
[]
```

**Key feature:** `block = '2 3'` ensures inertia is **never active** in the outer block.

#### 2. Radiation Damping - Region-Dependent

**Inner region (switchable):**
```
[Kernels]
  # Switched off during dynamic mode in inner region
  [rad_damp_x_inner]
    type = FarmsRadiationDamping
    variable = disp_x
    eta_constant = 1.85e7
    enable = true
    block = '2 3'
  []
  [rad_damp_y_inner]
    type = FarmsRadiationDamping
    variable = disp_y
    eta_constant = 1.85e7
    enable = true
    block = '2 3'
  []
[]
```

**Outer region (always active):**
```
[Kernels]
  # ALWAYS active - provides boundary wave absorption
  [rad_damp_x_outer]
    type = FarmsRadiationDamping
    variable = disp_x
    eta_constant = 1.85e7
    enable = true
    block = '1'
  []
  [rad_damp_y_outer]
    type = FarmsRadiationDamping
    variable = disp_y
    eta_constant = 1.85e7
    enable = true
    block = '1'
  []
[]
```

**Key feature:** Outer block radiation damping is **never disabled**, even during dynamic rupture.

#### 3. Sochacki Sponge (Outer Region)

Additional absorbing boundary treatment in block 1:

```
[Kernels]
  [sponge_damping_x]
    type = SochackiSpongeDamping
    variable = disp_x
    block = 1
    density = density
    sochacki_damping = sochacki_damping
    damping_scale = 2.0
  []
  [sponge_damping_y]
    type = SochackiSpongeDamping
    variable = disp_y
    block = 1
    density = density
    sochacki_damping = sochacki_damping
    damping_scale = 2.0
  []
[]

[Materials]
  [sochacki_sponge]
    type = SochackiSpongeMaterial
    block = 1
    inner_xmin = -20000
    inner_xmax = 20000
    inner_ymin = -20000
    inner_ymax = 20000
    outer_xmin = -30000
    outer_xmax = 30000
    outer_ymin = -30000
    outer_ymax = 30000
    s_max = 5.0
    profile = 'gaussian'
    gaussian_rate = 6.0
  []
[]
```

**Purpose:** Gradually increases damping from inner boundary to outer boundary, smoothly absorbing outgoing waves.

### Control System

```
[Controls]
  [strain_rate_switch]
    type = FarmsConditionalPostprocessorEnableControl
    postprocessor = max_dev_strain_rate
    comparison_type = greater_than
    threshold = 1e-5
    reverse_on_false = true

    # Enable when strain rate > 1e-5 (DYNAMIC mode in INNER region)
    enable_objects = 'Kernels/inertia_x Kernels/inertia_y
                      AuxKernels/accel_x_aux AuxKernels/accel_y_aux
                      AuxKernels/vel_x_aux_dynamic AuxKernels/vel_y_aux_dynamic'

    # Disable when strain rate > 1e-5 (turn off QUASI-DYNAMIC in INNER region only)
    disable_objects = 'Kernels/rad_damp_x_inner Kernels/rad_damp_y_inner
                       AuxKernels/vel_x_aux_quasi AuxKernels/vel_y_aux_quasi'
  []
[]
```

**Note:** Control only switches **inner region** kernels. Outer region remains untouched.

---

## Physical Behavior

### Outer Block (Block 1)

**Always quasi-dynamic:**
- Radiation damping equation: ∇·σ + η·v = 0
- No inertial effects (ρ·a = 0)
- Damping coefficient: η = 1.85e7 = 2μ/c_s
- Additional Sochacki sponge near boundaries

**Response to waves:**
1. Rupture waves propagate from inner region
2. Radiation damping dissipates wave energy
3. Sochacki sponge provides additional absorption
4. Waves decay before reflecting from boundaries

### Inner Region (Blocks 2, 3)

**Quasi-dynamic phase (ε̇ ≤ 1e-5):**
- Radiation damping active
- Loading accumulates stress
- No wave propagation
- dt can grow to 50 s

**Dynamic phase (ε̇ > 1e-5):**
- Full wave equation: ρ·a + ∇·σ = 0
- Inertia captures rupture dynamics
- Wave propagation resolved
- dt limited to 0.01 s

**Transition at boundary:**
- Smooth continuity at block interface
- Radiation damping in block 1 absorbs outgoing waves
- No spurious reflections from formulation change

---

## Advantages

### 1. Boundary Wave Isolation

**Problem prevented:**
- Loading BCs generate waves
- These waves would propagate with inertia in full dynamic
- Contaminate interior rupture dynamics

**Solution:**
- Outer block has no inertia
- Boundary waves are quasi-static disturbances
- Radiation damping dissipates them
- Clean rupture dynamics in inner region

### 2. Computational Efficiency

**Outer block:**
- Always quasi-dynamic → can use large dt during loading
- No need for small dt to resolve boundary waves
- Sponge provides additional stability

**Inner region:**
- Large dt (50 s) during loading
- Small dt (0.01 s) only during rupture
- Efficient overall time marching

### 3. Physical Consistency

**Outer region:**
- Represents far-field quasi-static response
- Physically reasonable for tectonic loading
- No artificial wave propagation from boundaries

**Inner region:**
- Full dynamic rupture when needed
- Natural transition to quasi-dynamic between events
- Captures both slow loading and fast slip

---

## Comparison with Alternatives

### Alternative 1: Full Dynamic Everywhere

```
Inertia active in all blocks at all times
```

**Problems:**
- Boundary loading generates inertial waves
- These waves propagate into domain
- Contaminate rupture dynamics
- Require very small dt even during loading
- Computationally expensive

### Alternative 2: Dashpot Absorbing Boundaries

```
Apply velocity-proportional damping at boundaries
```

**Issues:**
- Works for outgoing waves from rupture
- Doesn't prevent incoming waves from loading
- Still need to handle boundary condition application
- Can have impedance mismatch issues

### Alternative 3: Sponge Layer Only

```
Use only Sochacki sponge without inertia restriction
```

**Limitations:**
- Sponge absorbs waves but doesn't prevent their generation
- Boundary waves still carry inertia effects initially
- May need very thick sponge layer
- Less efficient than preventing inertia entirely

### Our Approach: Block-Selective Inertia + Sponge

```
Outer block: No inertia + radiation damping + sponge
Inner region: Adaptive inertia switching
```

**Benefits:**
- ✅ Prevents boundary wave generation (no inertia in outer block)
- ✅ Efficient loading phase (quasi-dynamic with large dt)
- ✅ Accurate rupture dynamics (full dynamic in inner region)
- ✅ Clean transition (radiation damping at interfaces)
- ✅ Additional wave absorption (Sochacki sponge)

---

## Verification

### Check 1: Block Identification

Verify mesh blocks are correctly assigned:
```bash
# In Paraview or from mesh file
- Block 1: Outer region (x: [-30k, -20k] ∪ [20k, 30k] or y: similar)
- Blocks 2, 3: Inner region (within [-20k, 20k] × [-20k, 20k])
```

### Check 2: Kernel Activity

Monitor which kernels are active during simulation:

**Loading phase:**
- `rad_damp_x_inner`, `rad_damp_y_inner`: ON (blocks 2, 3)
- `rad_damp_x_outer`, `rad_damp_y_outer`: ON (block 1)
- `inertia_x`, `inertia_y`: OFF

**Dynamic phase:**
- `inertia_x`, `inertia_y`: ON (blocks 2, 3 only)
- `rad_damp_x_inner`, `rad_damp_y_inner`: OFF (blocks 2, 3)
- `rad_damp_x_outer`, `rad_damp_y_outer`: STILL ON (block 1) ← **Critical**

### Check 3: Velocity Fields

Plot velocity magnitude over space and time:

**Expected:**
- **Outer block:** Very low velocities during loading (quasi-static)
- **Inner region (loading):** Low velocities (quasi-dynamic)
- **Inner region (rupture):** High velocities (dynamic)
- **Outer block (during rupture):** Moderate velocities from radiation damping, not inertial propagation

```python
import paraview.simple as pv

# Load exodus file
data = pv.OpenDataFile('unified_solve_main_eps1em4_cd10_cycle_out.e')

# Extract velocity magnitude
calc = pv.Calculator(Input=data)
calc.ResultArrayName = 'vel_mag'
calc.Function = 'sqrt(vel_x^2 + vel_y^2)'

# Clip to outer block
clip = pv.Clip(Input=calc)
clip.ClipType = 'Box'
clip.ClipType.Position = [-30000, -30000, 0]
clip.ClipType.Length = [60000, 60000, 0]
clip.Invert = 0

# Check maximum velocity in outer block during rupture
# Should be much lower than inner region peak velocities
```

### Check 4: Strain Rate Distribution

Verify strain rate criterion works correctly:

```bash
# Check CSV output
grep max_dev_strain_rate unified_solve_main_eps1em4_cd10_cycle_csv.csv

# Expected:
# - Loading: max_dev_strain_rate ~ 1e-7 to 1e-6 (below threshold)
# - Rupture: max_dev_strain_rate > 1e-5 (above threshold)
```

**Important:** Strain rate is monitored on **blocks 1 and 3** only (see AuxKernels), not the elastic outer block.

---

## Tuning Parameters

### Radiation Damping Coefficient

**Current:** η = 1.85e7 = 2μ/c_s

**To adjust:**
```
η = 2 × μ / c_s
  = 2 × 32.04e9 / 3464
  = 1.85e7 Pa·s/m
```

**Effect:**
- Higher η → stronger damping → waves dissipate faster
- Lower η → weaker damping → more wave propagation
- Typical: η = 2μ/c_s for impedance matching

### Sochacki Sponge Parameters

**Current:**
```
s_max = 5.0              # Maximum attenuation
sponge_profile = 'gaussian'
gaussian_rate = 6.0      # Steepness
```

**Tuning:**
- `s_max`: Increase to 10-20 for stronger absorption
- `gaussian_rate`: Increase for more rapid ramp-up
- `profile`: Options: 'linear', 'exponent', 'cubic', 'exponential', 'gaussian'

### Block Interface Location

**Current:**
- Inner region: [-20 km, 20 km] × [-20 km, 20 km]
- Outer region: [-30 km, 30 km] × [-30 km, 30 km]
- Buffer width: 10 km

**To adjust:** Modify mesh geometry in `mesh_test_outerblock.geo`
```
small_xmin = -20000  # Inner boundary
big_xmin = -30000    # Outer boundary
→ Buffer width = 10 km
```

**Guidelines:**
- Minimum buffer: ~3-5 wavelengths of expected rupture waves
- Wavelength: λ ~ c_s / f ~ 3464 / 1 Hz ~ 3.5 km
- Current 10 km buffer ~ 3 wavelengths (adequate)
- Increase to 15-20 km for extra safety

---

## Integration with Other Features

### Works with Adaptive Time Stepping

Block-selective inertia is fully compatible with:
- `FarmsIterationAdaptiveDT` timestepper
- `FarmsAdaptiveTimeStepBound` for mode-dependent dt_max
- Quasi-dynamic: dt → 50 s (in both inner and outer regions)
- Dynamic: dt → 0.01 s (controlled by inner region dynamics)

### Works with Damage-Breakage Model

Inner region uses damage-breakage constitutive model:
- Blocks 1, 3: `ComputeLagrangianDamageBreakageStressPK2Diffused`
- Block 2: Simple elastic (`ComputeStVenantKirchhoffStress`)

Outer block (block 1) can have damage if needed, or remain elastic for efficiency.

### Works with Multi-App System

Sub-app for damage evolution continues to work:
- Transfers occur for all blocks
- Damage evolution restricted to blocks where stress material is defined

---

## Common Issues

### Issue: Waves still appear in outer block

**Possible causes:**
1. Mesh blocks incorrectly assigned
2. Inertia kernels not properly restricted

**Debug:**
```bash
# Check kernel blocks in input file
grep -A 5 "inertia_x" unified_solve_main_eps1em4_cd10_cycle.i
# Should see: block = '2 3'
```

### Issue: Interface artifacts at block boundary

**Possible causes:**
1. Radiation damping coefficient mismatch
2. Material property discontinuity

**Solutions:**
1. Ensure consistent η in all blocks
2. Use continuous mesh across interface (no gaps)
3. Check that stress computation is consistent

### Issue: Simulation slower than expected

**Possible causes:**
1. Time step not growing during loading
2. Outer block still triggering small dt

**Check:**
```bash
# Monitor time step evolution
grep "_dt," unified_solve_main_eps1em4_cd10_cycle_csv.csv | tail -100

# Should see dt → 50 s during quasi-dynamic loading
# Should see dt → 0.01 s during dynamic rupture
```

---

## Summary

This block-selective approach provides:

1. **Wave Isolation**: Outer block prevents boundary loading waves from contaminating interior
2. **Efficiency**: Large time steps during loading (quasi-dynamic everywhere)
3. **Accuracy**: Full dynamic rupture in fault region when needed
4. **Stability**: Radiation damping + sponge absorb outgoing waves
5. **Flexibility**: Easy to adjust buffer zone size and damping parameters

The key insight: **Not everywhere needs inertia**. By restricting inertia to the region where rupture occurs, we get clean dynamics without boundary artifacts, while maintaining computational efficiency.

---

## Related Documentation

- `SWITCHING_README.md` - Dynamic/quasi-dynamic mode switching
- `ADAPTIVE_TIMESTEPPING_README.md` - Time step adaptation
- `unified_solve_main_eps1em4_cd10_cycle.i` - Full input file
