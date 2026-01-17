# Implementation Guide: Three-Field AD Formulation with Damage-Dependent Hydraulic Properties

## Objective

Compare the difference between two-field and three-field formulation on pulse power fracturing:
1. AD (Automatic Differentiation) compatibility
2. Phase-field damage coupling
3. Damage-dependent hydraulic properties evolution

## Current State Analysis

| Feature | Three-Field (existing) | Two-Field (reference) | Pure Solid AD (reference) |
|---------|------------------------|----------------------|---------------------------|
| Dimension | 3D | 2D | 2D |
| AD Support | Partial (AD kernels, but smeared cracking) | Non-AD | Full AD |
| Damage Model | Smeared cracking | Phase-field | Phase-field |
| Hydraulic Evolution | Fixed permeability from smeared cracking | Full damage-dependent (α, φ, M, k) | N/A |
| Variables | (disp_x,y,z), (wf_x,y,z), p | (disp_x,y), pp | (disp_x,y) |

### Reference Files

- **Three-field (existing)**: `pulsepower/pf_code2d_porousflow/threefield/elasticity.i`
- **Two-field (reference)**: `pulsepower/pf_code2d_porousflow/parametric_study/case_cf1_domain1x/elasticity.i`
- **Pure solid AD (reference)**: `pulsepower/pf_code2d_puresolid/code_engbudget/elasticity.i`

---

## Phase 1: Core Infrastructure

### Task 1.1: Create AD Phase-Field Elasticity Model for Three-Field

**Files to create**:
- `include/materials/phasefield_smalldeform_AD/ADSmallDeformationIsotropicElasticityThreeField.h`
- `src/materials/phasefield_smalldeform_AD/ADSmallDeformationIsotropicElasticityThreeField.C`

**Purpose**: Extend the existing `NDSmallDeformationIsotropicElasticity` pattern to AD version with three-field compatibility.

**Key features**:
- AD-compatible stress computation with spectral decomposition
- Damage-dependent permeability (Darcy-Poiseuille model)
- Output: `psie_active`, `effective_perm`, `solid_bulk_compliance_damaged`

**Reference**: `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C`

---

### Task 1.2: Create AD Damage-Dependent Hydraulic Properties Materials

#### Task 1.2a: `ElkADPorousFlowDamagedBiotCoefficient`

**Files to create**:
- `include/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotCoefficient.h`
- `src/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotCoefficient.C`

**Formula**: `α(d) = 1 - g(d)·K₀/Kₛ` where `g(d) = max((1-d)², η)`

**Output property**: `biot_coefficient_damaged`

**Reference**: `src/materials/porousflowmatprops/ElkPorousFlowDamagedBiotCoefficient.C`

---

#### Task 1.2b: `ElkADPorousFlowDamagedPorosity`

**Files to create**:
- `include/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.h`
- `src/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.C`

**Formula**: `φ(d) = φ₀ + (1-φ₀)·(1-g(d))`

**Output property**: `PorousFlow_porosity_qp_damaged`

**Reference**: `src/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.C`

---

#### Task 1.2c: `ElkADPorousFlowDamagedBiotModulus`

**Files to create**:
- `include/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotModulus.h`
- `src/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotModulus.C`

**Formula**: `1/M = φ(d)/Kf + (α(d)-φ(d))/Kₛ`

**Output property**: `biot_modulus_damaged`

**Reference**: `src/materials/porousflowmatprops/ElkPorousFlowDamagedBiotModulus.C`

---

### Task 1.3: Create `ElkADPoroDynamicPhaseFieldMaterials`

**Files to create**:
- `include/materials/threefieldporodynamics/ElkADPoroDynamicPhaseFieldMaterials.h`
- `src/materials/threefieldporodynamics/ElkADPoroDynamicPhaseFieldMaterials.C`

**Purpose**: Assemble all damage-dependent properties for three-field formulation.

**Key modifications from `ElkADPoroDynamicSmearedCrackingMaterials`**:
```cpp
// Current (fixed values):
_biot_coefficient[_qp] = _biot_coefficient_val;
_porosity[_qp] = _porosity_val;

// New (damage-dependent):
_biot_coefficient[_qp] = _biot_coefficient_damaged[_qp];
_porosity[_qp] = _porosity_damaged[_qp];
_biot_modulus[_qp] = _biot_modulus_damaged[_qp];
```

---

## Phase 2: Kernel Modifications

### Task 2.1: Create `ElkADDynamicDarcyFlowPhaseField`

**Files to create**:
- `include/kernels/threefieldporodynamics/ElkADDynamicDarcyFlowPhaseField.h`
- `src/kernels/threefieldporodynamics/ElkADDynamicDarcyFlowPhaseField.C`

**Purpose**: Use damaged permeability/porosity from material properties.

**Reference**: `ElkADDynamicDarcyFlowSmearedCracking`

---

### Task 2.2: Modify `ElkADPoroMechanicsCoupling`

Add option `use_damaged_biot = true` to use damaged Biot coefficient.

---

### Task 2.3: Modify `ElkADMassConservationNewmark`

Add option `use_damaged_biot = true` to use damaged Biot coefficient.

---

## Phase 3: Input File & Integration

### Task 3.1: Create 2D Input File Template

**File**: `pulsepower/pf_threefield/elasticity.i`

Structure with phase-field coupling and damage-dependent hydraulic properties.

---

## Phase 4: Validation

### Task 4.1: Energy Budget Implementation

Track energy components for comparison with two-field formulation.

---

## Unit Tests

### Phase 1 Unit Tests

#### Test 1.1: AD Damaged Biot Coefficient

**File**: `tests/phase1_materials/test_ad_damaged_biot_coefficient.i`

**Test Description**: Verify that the damaged Biot coefficient is computed correctly for different damage values.

**Test Cases**:
| Damage (d) | g(d) = (1-d)² | K₀ = 50e9 | Kₛ = 50e9 | α = 1 - g·K₀/Kₛ |
|------------|---------------|-----------|-----------|-----------------|
| 0.0        | 1.0           | 50e9      | 50e9      | 0.0             |
| 0.5        | 0.25          | 50e9      | 50e9      | 0.75            |
| 1.0        | η ≈ 0         | 50e9      | 50e9      | ~1.0            |

**Expected Behavior**:
- At d=0: α = 0 (intact material)
- At d=1: α → 1 (fully damaged, approaching incompressible)
- α should be clamped to [0, 1]

**Test Spec** (`tests/phase1_materials/tests`):
```
[Tests]
  [test_ad_damaged_biot_coefficient]
    type = 'CSVDiff'
    input = 'test_ad_damaged_biot_coefficient.i'
    csvdiff = 'test_ad_damaged_biot_coefficient_out.csv'
    requirement = 'AD damaged Biot coefficient shall be computed as α(d) = 1 - g(d)·K₀/Kₛ'
  []
  [test_biot_coefficient_bounds]
    type = 'RunApp'
    input = 'test_ad_damaged_biot_coefficient.i'
    cli_args = 'AuxVariables/d/initial_condition=1.5'
    expect_out = 'biot_coefficient_damaged'
    requirement = 'AD damaged Biot coefficient shall be clamped to [0, 1]'
  []
[]
```

---

#### Test 1.2: AD Damaged Porosity

**File**: `tests/phase1_materials/test_ad_damaged_porosity.i`

**Test Description**: Verify that the damaged porosity is computed correctly.

**Test Cases**:
| Damage (d) | g(d) = (1-d)² | φ₀ = 0.008 | φ(d) = φ₀ + (1-φ₀)·(1-g) |
|------------|---------------|------------|---------------------------|
| 0.0        | 1.0           | 0.008      | 0.008                     |
| 0.5        | 0.25          | 0.008      | 0.752                     |
| 1.0        | 0.0           | 0.008      | 1.0 (clamped to upper)    |

**Expected Behavior**:
- At d=0: φ = φ₀ (initial porosity)
- At d=1: φ → 1.0 (fully damaged becomes fully porous)
- φ should be clamped to [porosity_lower_bound, porosity_upper_bound]

**Test Spec**:
```
[Tests]
  [test_ad_damaged_porosity]
    type = 'CSVDiff'
    input = 'test_ad_damaged_porosity.i'
    csvdiff = 'test_ad_damaged_porosity_out.csv'
    requirement = 'AD damaged porosity shall be computed as φ(d) = φ₀ + (1-φ₀)·(1-g(d))'
  []
[]
```

---

#### Test 1.3: AD Damaged Biot Modulus

**File**: `tests/phase1_materials/test_ad_damaged_biot_modulus.i`

**Test Description**: Verify that the damaged Biot modulus is computed correctly.

**Formula**: `1/M = φ(d)/Kf + (α(d)-φ(d))/Kₛ`

**Test Cases**:
| d   | α(d)  | φ(d)  | Kf = 2.24e9 | Kₛ = 50e9 | 1/M                  | M        |
|-----|-------|-------|-------------|-----------|----------------------|----------|
| 0.0 | 0.0   | 0.008 | 2.24e9      | 50e9      | 3.57e-12 + (-0.008)/Kₛ | ~2.8e11  |
| 0.5 | 0.75  | 0.752 | 2.24e9      | 50e9      | 3.36e-10 + (-0.002)/Kₛ | ~2.98e9  |

**Test Spec**:
```
[Tests]
  [test_ad_damaged_biot_modulus]
    type = 'CSVDiff'
    input = 'test_ad_damaged_biot_modulus.i'
    csvdiff = 'test_ad_damaged_biot_modulus_out.csv'
    requirement = 'AD damaged Biot modulus shall be computed as 1/M = φ/Kf + (α-φ)/Kₛ'
  []
[]
```

---

#### Test 1.4: AD Elasticity with Permeability

**File**: `tests/phase1_materials/test_ad_elasticity_permeability.i`

**Test Description**: Verify that the effective permeability is computed correctly using Darcy-Poiseuille model.

**Formula**: `k(d) = k₀ + d^n · (w²/12 - k₀)` where `w = d·wc`

**Test Cases**:
| d   | k₀ = 5e-19 | wc = 2.9e-4 | n = 10 | k(d)         |
|-----|------------|-------------|--------|--------------|
| 0.0 | 5e-19      | 2.9e-4      | 10     | 5e-19        |
| 0.5 | 5e-19      | 2.9e-4      | 10     | ~5.87e-19    |
| 1.0 | 5e-19      | 2.9e-4      | 10     | ~7.0e-9      |

**Test Spec**:
```
[Tests]
  [test_ad_elasticity_permeability]
    type = 'CSVDiff'
    input = 'test_ad_elasticity_permeability.i'
    csvdiff = 'test_ad_elasticity_permeability_out.csv'
    requirement = 'Effective permeability shall follow Darcy-Poiseuille model'
  []
[]
```

---

#### Test 1.5: AD PoroDynamic Materials Assembly

**File**: `tests/phase1_materials/test_ad_porodynamic_materials.i`

**Test Description**: Verify that all materials are assembled correctly in `ElkADPoroDynamicPhaseFieldMaterials`.

**Test Cases**:
- Verify density = ρs·(1-φ) + ρf·φ
- Verify Biot modulus uses damaged α and φ
- Verify effective permeability is passed through

**Test Spec**:
```
[Tests]
  [test_ad_porodynamic_materials_density]
    type = 'CSVDiff'
    input = 'test_ad_porodynamic_materials.i'
    csvdiff = 'test_ad_porodynamic_materials_out.csv'
    requirement = 'Density shall be computed as ρ = ρs·(1-φ) + ρf·φ with damaged porosity'
  []
[]
```

---

### Phase 2 Unit Tests

#### Test 2.1: AD Dynamic Darcy Flow

**File**: `tests/phase2_kernels/test_ad_dynamic_darcy_flow.i`

**Test Description**: Verify that the dynamic Darcy flow kernel uses damaged properties correctly.

**Test Spec**:
```
[Tests]
  [test_ad_dynamic_darcy_flow]
    type = 'Exodiff'
    input = 'test_ad_dynamic_darcy_flow.i'
    exodiff = 'test_ad_dynamic_darcy_flow_out.e'
    requirement = 'AD dynamic Darcy flow kernel shall use damaged permeability'
  []
[]
```

---

#### Test 2.2: AD PoroMechanics Coupling with Damaged Biot

**File**: `tests/phase2_kernels/test_ad_poromech_coupling.i`

**Test Description**: Verify that the poromechanics coupling uses damaged Biot coefficient.

**Test Spec**:
```
[Tests]
  [test_ad_poromech_coupling]
    type = 'Exodiff'
    input = 'test_ad_poromech_coupling.i'
    exodiff = 'test_ad_poromech_coupling_out.e'
    requirement = 'AD poromechanics coupling kernel shall use damaged Biot coefficient when enabled'
  []
[]
```

---

### Phase 3 Integration Tests

#### Test 3.1: Full Three-Field with Phase-Field

**File**: `tests/phase3_integration/test_threefield_phasefield.i`

**Test Description**: Run a complete three-field simulation with phase-field damage.

**Verification**:
- Energy balance: input energy ≈ kinetic + elastic + fracture + dissipation
- Damage evolution follows phase-field model
- Hydraulic properties evolve with damage

**Test Spec**:
```
[Tests]
  [test_threefield_phasefield]
    type = 'CSVDiff'
    input = 'test_threefield_phasefield.i'
    csvdiff = 'test_threefield_phasefield_out.csv'
    rel_err = 1e-5
    requirement = 'Three-field phase-field simulation shall conserve energy'
  []
[]
```

---

#### Test 3.2: Two-Field vs Three-Field Comparison

**File**: `tests/phase3_integration/test_comparison.i`

**Test Description**: Compare results between two-field and three-field formulations.

**Metrics**:
- Damage pattern similarity
- Pressure field comparison
- Energy partitioning

**Test Spec**:
```
[Tests]
  [test_twofield_threefield_comparison]
    type = 'CSVDiff'
    input = 'test_comparison.i'
    csvdiff = 'test_comparison_out.csv'
    requirement = 'Two-field and three-field shall produce comparable damage patterns'
  []
[]
```

---

## Implementation Sequence

```
Phase 1: Core Infrastructure
├── 1.1 ADSmallDeformationIsotropicElasticityThreeField
│   └── Compute: stress, psie_active, effective_perm
│
├── 1.2 AD Damaged Properties Materials
│   ├── 1.2a ElkADPorousFlowDamagedBiotCoefficient
│   ├── 1.2b ElkADPorousFlowDamagedPorosity
│   └── 1.2c ElkADPorousFlowDamagedBiotModulus
│
└── 1.3 ElkADPoroDynamicPhaseFieldMaterials
    └── Assemble all damaged properties

Phase 2: Kernel Modifications
├── 2.1 ElkADDynamicDarcyFlowPhaseField
├── 2.2 Modify ElkADPoroMechanicsCoupling
└── 2.3 Modify ElkADMassConservationNewmark

Phase 3: Input File & Integration
├── 3.1 Create 2D input file template
├── 3.2 Setup MultiApp for phase-field
└── 3.3 Configure mesh

Phase 4: Validation
├── 4.1 Energy budget tracking
├── 4.2 Compare with two-field results
└── 4.3 Verify damage-hydraulic coupling
```

---

## Key Property Dependencies

```
                    +---------------+
                    | Phase-field   |
                    |     d         |
                    +-------+-------+
                            |
            +---------------+---------------+
            |               |               |
            v               v               v
     +-----------+   +-----------+   +-----------+
     | g(d)      |   | effective |   | solid_bulk|
     | =(1-d)^2  |   | permeab.  |   | compliance|
     +-----+-----+   +-----+-----+   +-----+-----+
           |               |               |
     +-----+-----+         |         +-----+-----+
     |           |         |         |           |
     v           v         v         v           v
+--------+  +--------+ +--------+ +--------+ +--------+
| alpha  |  | phi    | | k(d)   | | K(d)   | | M(d)   |
| Biot   |  | poros- | | perme- | | bulk   | | Biot   |
| coeff  |  | ity    | | ability| | modulus| | modulus|
+----+---+  +----+---+ +----+---+ +--------+ +----+---+
     |           |          |                      |
     +-----------+----------+----------------------+
                            |
                            v
               +------------------------+
               | ElkADPoroDynamic       |
               | PhaseFieldMaterials    |
               +------------------------+
```

---

## Summary of Tasks

| # | Task | Priority | Complexity | Dependencies |
|---|------|----------|------------|--------------|
| 1.1 | AD Phase-Field Elasticity Model | High | Medium | None |
| 1.2a | AD Damaged Biot Coefficient | High | Low | 1.1 |
| 1.2b | AD Damaged Porosity | High | Low | 1.1 |
| 1.2c | AD Damaged Biot Modulus | High | Low | 1.2a, 1.2b |
| 1.3 | AD PoroDynamic PhaseField Materials | High | Medium | 1.2a-c |
| 2.1 | AD Dynamic Darcy Flow PhaseField | High | Medium | 1.3 |
| 2.2 | Modify PoroMechanicsCoupling | Medium | Low | 1.2a |
| 2.3 | Modify MassConservationNewmark | Medium | Low | 1.2a |
| 3.1 | Create 2D input file template | High | Low | All above |
| 4.1 | Energy budget implementation | Medium | Low | 3.1 |

---

## Progress Status

### Phase 1: Core Infrastructure - **COMPLETED**

#### Completed Tasks:

**Task 1.2a: ElkADPorousFlowDamagedBiotCoefficient** ✅
- Header: `include/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotCoefficient.h`
- Source: `src/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotCoefficient.C`
- Unit test: `tests/phase1_materials/test_ad_damaged_biot_coefficient.i`

**Task 1.2b: ElkADPorousFlowDamagedPorosity** ✅
- Header: `include/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.h`
- Source: `src/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.C`
- Unit test: `tests/phase1_materials/test_ad_damaged_porosity.i`

**Task 1.2c: ElkADPorousFlowDamagedBiotModulus** ✅
- Header: `include/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotModulus.h`
- Source: `src/materials/porousflowmatprops/ElkADPorousFlowDamagedBiotModulus.C`
- Unit test: `tests/phase1_materials/test_ad_damaged_biot_modulus.i`

**Task 1.3: ElkADPoroDynamicPhaseFieldMaterials** ✅
- Header: `include/materials/threefieldporodynamics/ElkADPoroDynamicPhaseFieldMaterials.h`
- Source: `src/materials/threefieldporodynamics/ElkADPoroDynamicPhaseFieldMaterials.C`
- Unit test: `tests/phase1_materials/test_ad_porodynamic_materials.i`

**Test specification file** ✅
- `tests/phase1_materials/tests` - Contains all four Phase 1 unit tests

### To Run All Phase 1 Tests:

```bash
cd /path/to/moose
./run_tests -i pulsepower/pf_threefield/tests/phase1_materials
```

Or run individual tests:
```bash
./run_tests -i pulsepower/pf_threefield/tests/phase1_materials/test_ad_damaged_biot_coefficient.i
./run_tests -i pulsepower/pf_threefield/tests/phase1_materials/test_ad_damaged_porosity.i
./run_tests -i pulsepower/pf_threefield/tests/phase1_materials/test_ad_damaged_biot_modulus.i
./run_tests -i pulsepower/pf_threefield/tests/phase1_materials/test_ad_porodynamic_materials.i
```

### Phase 2: Kernel Modifications - **COMPLETED**

#### Completed Tasks:

**Task 2.1: ElkADDynamicDarcyFlowPhaseField** ✅
- Header: `include/kernels/threefieldporodynamics/ElkADDynamicDarcyFlowPhaseField.h`
- Source: `src/kernels/threefieldporodynamics/ElkADDynamicDarcyFlowPhaseField.C`
- Unit test: `tests/phase2_kernels/test_ad_dynamic_darcy_flow_phasefield.i`
- **Key change**: Uses `RankTwoTensor` for permeability (AD-compatible with phase-field materials)

**Task 2.2: ElkADPoroMechanicsCoupling** ✅ (No modification needed)
- The kernel already gets `biot_coefficient` from material properties
- When `multiply_biot_coefficient = true`, it uses the damaged Biot coefficient
- Unit test: `tests/phase2_kernels/test_ad_poromechanics_coupling_damaged.i`

**Task 2.3: ElkADMassConservationNewmark** ✅ (No modification needed)
- The kernel already gets `biot_coefficient` and `biot_modulus` from material properties
- When `multiply_biot_coefficient = true`, it uses the damaged properties
- Unit test: `tests/phase2_kernels/test_ad_mass_conservation_damaged.i`

**Test specification file** ✅
- `tests/phase2_kernels/tests` - Contains all three Phase 2 unit tests

### To Run All Phase 2 Tests:

```bash
cd /path/to/moose
./run_tests --spec-file pulsepower/pf_threefield/tests/phase2_kernels/tests
```

Or run individual tests:
```bash
./run_tests -i pulsepower/pf_threefield/tests/phase2_kernels/test_ad_dynamic_darcy_flow_phasefield.i
./run_tests -i pulsepower/pf_threefield/tests/phase2_kernels/test_ad_poromechanics_coupling_damaged.i
./run_tests -i pulsepower/pf_threefield/tests/phase2_kernels/test_ad_mass_conservation_damaged.i
```

### Phase 3: Input File & Integration - **COMPLETED**

#### Completed Tasks:

**Task 3.1: 2D Input File Template** ✅
- File: `pulsepower/pf_threefield/elasticity_phasefield.i`
- Full three-field (u, w, p) formulation with phase-field damage
- MultiApp coupling to `fracture.i` for phase-field evolution
- AD damage-dependent hydraulic properties
- Energy tracking postprocessors

**Task 3.2: Integration Test** ✅
- File: `tests/phase3_integration/test_threefield_phasefield_2d.i`
- Simplified test with prescribed damage (no MultiApp)
- Verifies three-field equations work with damage-dependent properties

**Test specification file** ✅
- `tests/phase3_integration/tests`

### To Run Phase 3 Tests:

```bash
./run_tests --spec-file pulsepower/pf_threefield/tests/phase3_integration/tests
```

### Key Differences: Two-Field vs Three-Field

| Aspect | Two-Field (u, p) | Three-Field (u, w, p) |
|--------|------------------|----------------------|
| Variables | disp_x, disp_y, pp | disp_x, disp_y, wf_x, wf_y, p |
| Fluid Kinetic Energy | Derived from Darcy velocity | Direct from fluid velocity vf |
| Darcy Flow Kernel | PorousFlowFullySaturatedDarcyBase | ElkADDynamicDarcyFlowPhaseField |
| Fluid Inertia | Implicit in mass balance | Explicit: ElkADPoreFluidInertialForceCoupling |

### Phase 4: Validation - **PENDING**

Next steps:
1. Run comparison between two-field and three-field on same problem
2. Verify energy conservation
3. Compare damage patterns and pressure distributions
