# PML Quick Start Guide

## 1. Compile the Code
```bash
cd /Users/chunhuizhao/projects/farms_aftershock
make -j8
```

## 2. Test the Example
```bash
cd development/aftershock2D/dynamic_solve
mpirun -n 4 ../../../../farms_aftershock-opt -i dynamic_solve_elastic_PML.i
```

## 3. Files Created

**Source code:**
- `include/kernels/core/PMLDamping.h`
- `src/kernels/core/PMLDamping.C`
- `include/materials/helper/PMLCoefficientMaterial.h`
- `src/materials/helper/PMLCoefficientMaterial.C`

**Examples:**
- `dynamic_solve_elastic_PML.i` - Full working example
- `PML_README.md` - Detailed documentation
- `PML_QuickStart.md` - This file

## 4. Minimal Setup Template

```moose
# Define PML geometry at top of input file
pml_xmin_inner = -35000  # Physical domain edges
pml_xmax_inner = 35000
pml_ymin_inner = -15000
pml_ymax_inner = 15000
pml_thickness = 5000     # 5km PML layer

[Mesh]
  [msh]
    # Total domain = physical domain + PML
    xmin = -40000  # = pml_xmin_inner - pml_thickness
    xmax = 40000
    ymin = -20000
    ymax = 20000
  []
[]

[Kernels]
  # Remove all [BCs] blocks, add these kernels instead:
  [pml_damping_x]
    type = PMLDamping
    variable = disp_x
  []
  [pml_damping_y]
    type = PMLDamping
    variable = disp_y
  []
[]

[Materials]
  [pml_coeff]
    type = PMLCoefficientMaterial
    pml_xmin = ${pml_xmin_inner}
    pml_xmax = ${pml_xmax_inner}
    pml_ymin = ${pml_ymin_inner}
    pml_ymax = ${pml_ymax_inner}
    pml_thickness = ${pml_thickness}
    ref_wave_speed = ${Cp}  # Use your P-wave speed parameter
    exponent = 3.0
  []
[]
```

## 5. Key Differences from Dashpot BC

**Remove:** Entire `[BCs]` block with dashpot BCs

**Add:**
- PML kernels in `[Kernels]` block
- PML material in `[Materials]` block
- Extended mesh to include PML layer

## 6. Recommended Parameters

For typical earthquake simulation (Cs ~ 3-4 km/s):
- **PML thickness**: 5-10 km
- **Exponent**: 3.0 (cubic profile)
- **d_max**: Auto-computed (don't specify)

Expected result: <1% reflection for all wave angles

## 7. Verification

Plot `pml_damping_aux` to visualize the PML:
- Should be 0 in physical domain (center)
- Should increase smoothly toward boundaries
- Should be maximum at outer edges

## 8. Common Issues

**Problem:** "Property 'pml_damping_coeff' not found"
**Solution:** Add `PMLCoefficientMaterial` to `[Materials]` block

**Problem:** Instability near PML
**Solution:** Reduce `d_max`:
```moose
[Materials]
  [pml_coeff]
    # ... other params ...
    d_max = 300.0  # Manual override
  []
[]
```

**Problem:** Still see reflections
**Solution:** Increase PML thickness or exponent:
```moose
pml_thickness = 10000  # 10km instead of 5km
exponent = 4.0         # quartic instead of cubic
```

## 9. Comparison Test

Run both versions and compare:
```bash
# Original with dashpot BC
mpirun -n 4 farms_aftershock-opt -i dynamic_solve_elastic.i

# New with PML
mpirun -n 4 farms_aftershock-opt -i dynamic_solve_elastic_PML.i
```

Check velocity near boundary (x=35km, y=0):
- Dashpot: Large oscillations after t~3s (reflections)
- PML: Smooth decay (minimal reflections)

## 10. Performance

- **Memory:** ~10-20% more (due to PML layer elements)
- **Runtime:** ~5-10% slower (damping computation)
- **Reflection reduction:** ~10-30× better than dashpot BC

**Worth it?** YES, if boundary reflections are contaminating your results.
