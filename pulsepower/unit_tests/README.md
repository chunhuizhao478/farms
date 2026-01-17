# Unit Tests for FARMS Phase Field Implementations

## Overview

This directory contains unit tests for:
1. Rate-dependent phase field fracture (Hofacker & Miehe 2012)
2. Three-field porodynamics with phase field damage

## Directory Structure

```
unit_tests/
├── kernels/                    # Rate-dependent phase field kernel tests
├── phase1_materials/           # Three-field material property tests
├── phase2_kernels/             # Three-field kernel tests
├── phase3_integration/         # Three-field integration tests
├── run_tests.sh               # Test runner script
└── README.md                  # This file
```

---

## 1. Rate-Dependent Phase Field Tests (`kernels/`)

Based on Hofacker & Miehe (2012) - IJNME 93:276-301

### Test Files

| File | Description |
|------|-------------|
| `test_ADPFFViscousResistance.i` | Basic functionality (eta = 1e-6) |
| `test_ADPFFViscousResistance_zero_eta.i` | Rate-independent limit (eta = 0) |
| `test_ADPFFViscousResistance_high_eta.i` | High viscosity effect (eta = 1e-3) |

### Running Tests

```bash
# From farms_cdms directory
./farms-opt -i pulsepower/unit_tests/kernels/test_ADPFFViscousResistance.i
```

---

## 2. Three-Field Porodynamics Material Tests (`phase1_materials/`)

Tests for damaged porous media material properties.

### Test Files

| File | Description |
|------|-------------|
| `test_ad_damaged_biot_coefficient.i` | Tests damaged Biot coefficient: alpha_d = alpha_0 * g(d) |
| `test_ad_damaged_biot_modulus.i` | Tests damaged Biot modulus: M_d = M_0 * g(d) |
| `test_ad_damaged_porosity.i` | Tests damaged porosity: phi_d = phi_0 + (alpha-phi_0)*(1-g(d)) |
| `test_ad_porodynamic_materials.i` | Combined material property test |
| `test_ad_material_property_integration.i` | Material property integration test |

### Running Tests

```bash
# From farms_cdms directory
./farms-opt -i pulsepower/unit_tests/phase1_materials/test_ad_damaged_biot_coefficient.i
./farms-opt -i pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i
```

---

## 3. Three-Field Porodynamics Kernel Tests (`phase2_kernels/`)

Tests for coupled poromechanics and fluid flow kernels with damage.

### Test Files

| File | Description |
|------|-------------|
| `test_ad_dynamic_darcy_flow_phasefield.i` | Dynamic Darcy flow with phase field damage |
| `test_ad_mass_conservation_damaged.i` | Mass conservation in damaged porous media |
| `test_ad_poromechanics_coupling_damaged.i` | Poromechanics coupling with damage |

### Running Tests

```bash
# From farms_cdms directory
./farms-opt -i pulsepower/unit_tests/phase2_kernels/test_ad_dynamic_darcy_flow_phasefield.i
./farms-opt -i pulsepower/unit_tests/phase2_kernels/test_ad_mass_conservation_damaged.i
```

---

## 4. Three-Field Integration Tests (`phase3_integration/`)

Full integration tests for the three-field (displacement, pressure, phase field) model.

### Test Files

| File | Description |
|------|-------------|
| `test_threefield_phasefield_2d.i` | 2D three-field phase field simulation |

### Running Tests

```bash
# From farms_cdms directory
./farms-opt -i pulsepower/unit_tests/phase3_integration/test_threefield_phasefield_2d.i
```

---

## Running All Tests

Use the provided test runner script:

```bash
cd pulsepower/unit_tests
./run_tests.sh ../../farms-opt
```

Or run individual test directories:

```bash
# From farms_cdms directory
for f in pulsepower/unit_tests/phase1_materials/*.i; do
    ./farms-opt -i "$f"
done
```

---

## Verification Criteria

### Rate-Dependent Tests
- eta = 0 should give rate-independent (quasi-static) results
- Higher eta should slow down damage evolution
- All tests should converge without errors

### Three-Field Material Tests
- Material properties should degrade correctly with damage
- Biot coefficient: alpha_d -> alpha_0 as d -> 0, alpha_d -> 0 as d -> 1
- Porosity should increase with damage

### Three-Field Kernel Tests
- Mass conservation should be satisfied
- Poromechanics coupling should be stable
- Damage should affect fluid flow appropriately

---

## Troubleshooting

### Convergence Issues
1. Reduce time step size (dt)
2. Increase nl_max_its
3. Check material property values

### Missing Output Files
1. Ensure exodus = true in Outputs block
2. Check file permissions

---

## References

1. Hofacker, M., & Miehe, C. (2013). A phase field model of dynamic fracture. IJNME, 93(3), 276-301.
2. Miehe, C., et al. (2015). Phase field modeling of fracture in porous media. CMAME.
