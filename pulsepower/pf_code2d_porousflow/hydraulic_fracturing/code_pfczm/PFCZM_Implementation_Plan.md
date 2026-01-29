# PF-CZM Implementation Plan for Hydraulic Fracturing Simulation

## 1. Overview

This document outlines the implementation plan for transitioning from the current **AT1 phase-field model** to the **Phase-Field Cohesive Zone Model (PF-CZM)** based on Wu's unified phase-field theory (JMPS 2017).

### 1.1 Current Implementation (AT1)
- **Crack geometric function**: `α(d) = d` (linear)
- **Degradation function**: `g(d) = (1-d)^2*(1-η) + η` (quadratic power function)
- **Normalization constant**: `c0 = 8/3`
- **Characteristic**: No explicit failure strength; fracture initiates immediately under load

### 1.2 Target Implementation (PF-CZM)
- **Crack geometric function**: `α(d) = 2d - d²` (optimal for quasi-brittle materials)
- **Degradation function**: Rational function with explicit failure strength
- **Normalization constant**: `c0 = π`
- **Characteristic**: Well-defined failure strength and softening behavior; converges to cohesive zone model

---

## 2. Theoretical Background

### 2.1 Wu's Unified Phase-Field Theory

From Wu (JMPS 2017), the key governing equations are:

**Regularized crack surface functional:**
```
Ad(d) = ∫_B γ(d,∇d) dV
γ(d,∇d) = (1/c0) * [α(d)/l + l|∇d|²]
```

**Optimal crack geometric function (Eq. 3.18a in Wu):**
```
α(d) = ξd + (1-ξ)d²    with ξ = 2
     = 2d - d²
```

**Normalization constant:**
```
c0 = 4∫₀¹ √α(β) dβ = π    (for ξ = 2)
```

**Optimal energetic degradation function (Eq. 3.18b in Wu):**
```
ω(d) = (1-d)^p / [(1-d)^p + Q(d)]

where Q(d) = a1 * d * P(d)
      P(d) = 1 + a2*d * (1 + a3*d)
```

**Key parameters:**
```
a1 = (4/π) * (lch/l) = (4/π) * (E*Gc)/(ft²*l)    [Eq. 3.19a]
a2 = 2*(-2*k0*Gf/ft²)^(2/3) - (p + 1/2)          [Eq. 3.19b]
a3 = function of wc, ft, Gf (for p=2)             [Eq. 3.19c]
```

where:
- `lch = E*Gc/ft²`: Irwin's characteristic length
- `ft`: Tensile strength
- `Gc`: Fracture energy (critical energy release rate)
- `k0`: Initial slope of softening curve
- `wc`: Ultimate crack opening displacement
- `p`: Power exponent (typically p=2)

### 2.2 Softening Laws

| Softening Law | k0 | wc | p | a2 | a3 |
|---------------|----|----|---|----|----|
| Linear | -ft²/(2Gf) | 2Gf/ft | 2 | -0.5 | 0 |
| Exponential | -ft²/Gf | ∞ | 5/2 | 0.1748 | 0 |
| Hyperbolic | -2ft²/Gf | ∞ | 4 | 0.5379 | 0 |
| Cornelissen | -1.3546ft²/Gf | 5.1361Gf/ft | 2 | 1.3868 | 0.6567 |

### 2.3 Comparison: AT1 vs PF-CZM

| Feature | AT1 | PF-CZM |
|---------|-----|--------|
| α(d) | d | 2d - d² |
| c0 | 8/3 | π |
| Failure strength | Undefined | ft = √(E*Gc*ξ/(c0*l*a1)) |
| Half bandwidth | 2l | πl/2 |
| Support | Infinite | Finite (bounded) |
| Length scale independence | No | Yes (for mode-I) |
| Softening control | Limited | Full (via a2, a3) |

---

## 3. Implementation Details (COMPLETED)

### 3.1 Key Implementation Decision

**IMPORTANT**: The `NDSmallDeformationIsotropicElasticity` material in MOOSE/RACCOON computes the degradation function `g(d)` **internally** based on `model_type`. For PF-CZM, we must use:

```
model_type = PF_CZM
```

This requires providing material properties `a1`, `a2`, `a3`, and `p` which the code uses to compute `g`, `dg_dd`, and `d2g_dd2` internally (see `NDSmallDeformationIsotropicElasticity.C` lines 394-441).

### 3.2 Elasticity Application (elasticityhf_pfczm.i / elasticityhf_static_pfczm.i)

#### 3.2.1 PF-CZM Parameters

```
# =============================================================================
# PF-CZM Specific Parameters
# =============================================================================
# Tensile strength - KEY PARAMETER for PF-CZM
ft = 3e6            # Tensile strength [Pa] - from material testing

# Irwin's characteristic length
lch = '${fparse E * Gc_const / ft^2}'

# Normalization constant for alpha(d) = 2d - d^2
c0_val = 3.14159265359  # pi

# Initial slope xi of crack geometric function: xi = d(alpha)/dd at d=0
# For alpha(d) = 2d - d^2: xi = 2
xi_val = 2

# PF-CZM constitutive parameters
# Linear softening: p=2, a2=-0.5, a3=0
p_deg = 2
a2 = -0.5
a3 = 0
eta = 1e-6

# Critical fracture energy for degradation function
psic = '${fparse 3.0 * Gc_const / (8.0 * l)}'

# a1 parameter
a1_check = '${fparse 4.0 / c0_val * E * Gc_const / (ft^2 * l)}'
```

#### 3.2.2 Material Properties for PF-CZM Parameters

```
[Materials]
  # PF-CZM Parameters as Material Properties
  # Note: NDSmallDeformationIsotropicElasticity computes g internally for PF_CZM
  [pfczm_params]
    type = ADGenericConstantMaterial
    prop_names = 'a1_mat a2_mat a3_mat p_mat'
    prop_values = '${a1_check} ${a2} ${a3} ${p_deg}'
  []

  # PF-CZM Material Properties for fracture sub-app
  [pfczm_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc psic xi c0'
    prop_values = '${l} ${Gc_const} ${psic} ${xi_val} ${c0_val}'
  []
[]
```

#### 3.2.3 Crack Geometric Function

```
[crack_geometric]
  type = CrackGeometricFunction
  property_name = alpha
  expression = '2*d - d^2'
  phase_field = d
[]
```

#### 3.2.4 Elasticity Model with PF-CZM

```
[elasticity]
  type = NDSmallDeformationIsotropicElasticity
  block = domain
  bulk_modulus = K
  shear_modulus = G
  phase_field = d
  strain_energy_density = psie
  strain_energy_density_active = psie_active
  strain_energy_density_inactive = psie_inactive
  strain_energy_density_derivative = dpsie_dd
  degradation_function = g
  degradation_function_derivative = dg_dd
  degradation_function_second_derivative = d2g_dd2
  decomposition = SPECTRAL

  # PF-CZM MODEL TYPE - uses internal rational degradation function
  model_type = PF_CZM
  a1 = a1_mat
  a2 = a2_mat
  a3 = a3_mat
  p = p_mat

  eta = ${eta}
  output_properties = 'elastic_strain psie_active'
  outputs = exodus

  # Porous flow coupling
  porous_flow_coupling = true
  darcy_poiseuille_permeability_model = true
  intrinsic_permeability = ${intrinsic_permeability_domain}
  wc = ${wc}
  perm_exponent = ${perm_exponent}
[]
```

### 3.3 Fracture Sub-Application (fracturehf_pfczm.i / fracturehf_static_pfczm.i)

#### 3.3.1 Material Properties

```
[fracture_properties]
  type = ADGenericConstantMaterial
  prop_names = 'l Gc psic xi c0'
  prop_values = '${l} ${Gc_const} ${psic} 2 3.14159265359'
[]
```

#### 3.3.2 Crack Geometric Function

```
[crack_geometric]
  type = CrackGeometricFunction
  property_name = alpha
  expression = '2*d - d^2'
  phase_field = d
[]
```

#### 3.3.3 Degradation Function (RationalDegradationFunction)

```
[degradation]
  type = RationalDegradationFunction
  property_name = g
  expression = (1-d)^p/((1-d)^p+(Gc/psic*xi/c0/l)*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
  phase_field = d
  material_property_names = 'Gc psic xi c0 l'
  parameter_names = 'p a2 a3 eta'
  parameter_values = '${p_deg} ${a2} ${a3} ${eta}'
[]
```

#### 3.3.4 Free Energy Functional

```
[psi]
  type = ADDerivativeParsedMaterial
  property_name = psi
  expression = 'alpha*Gc/c0/l+0.5*g*psie_active'
  coupled_variables = 'd psie_active'
  material_property_names = 'alpha(d) g(d) Gc c0 l'
  derivative_order = 1
[]
```

**Note**: The `0.5` factor is needed because `psie_active` already contains a factor of 2 from the strain energy computation.

---

## 4. Implementation Checklist

### Phase 1: Basic PF-CZM Implementation
- [x] Update crack geometric function: `α(d) = 2d - d²`
- [x] Implement RationalDegradationFunction for degradation (fracture sub-app)
- [x] Use `model_type = PF_CZM` in NDSmallDeformationIsotropicElasticity (elasticity app)
- [x] Update normalization constant c0 = π
- [x] Add psic (critical fracture energy) material property
- [x] Add a1, a2, a3, p as material properties for PF_CZM model
- [x] Update free energy expression with 0.5 factor

### Phase 2: Parameter Calibration
- [x] Define tensile strength ft parameter
- [x] Calculate derived parameters (a1, a2, a3) for linear softening
- [x] Verify psic computation: `psic = 3*Gc/(8*l)`
- [ ] Test with 1D uniaxial tension case

### Phase 3: Integration with Porous Flow
- [x] Verify compatibility with enhanced history energy
- [x] Check permeability model interaction (Darcy-Poiseuille)
- [x] Update Biot coefficient damage coupling
- [ ] Validate energy conservation

### Phase 4: Validation
- [ ] Compare with mode-I fracture benchmark (3-point bending)
- [ ] Verify length-scale independence
- [ ] Test softening curve recovery
- [ ] Run hydraulic fracturing simulations

---

## 5. Key Differences Summary

| Component | AT1 | PF-CZM |
|-----------|-----|--------|
| `alpha` expression | `d` | `2*d - d^2` |
| `c0` | 8/3 ≈ 2.667 | π ≈ 3.14159 |
| `model_type` in elasticity | `AT1` | `PF_CZM` |
| Degradation type (fracture) | PowerDegradationFunction | RationalDegradationFunction |
| `g` expression | `(1-d)^p*(1-η)+η` | `(1-d)^p/((1-d)^p+a1*d*P(d))*(1-η)+η` |
| New parameters | None | `ft`, `psic`, `a1`, `a2`, `a3`, `p` |
| Half bandwidth Du | 2l | πl/2 ≈ 1.57l |
| Failure strength | Undefined | ft = √(E·Gc·ξ/(c0·l·a1)) |

---

## 6. Source Code Implementation

### 6.1 NDSmallDeformationIsotropicElasticity.C

The PF-CZM degradation function is computed internally in `NDSmallDeformationIsotropicElasticity::computeGDerivatives()` (lines 394-441):

```cpp
else if (_model_type == "PF_CZM"){
  // Get the parameters
  const Real a1 = (*_a1_prop)[_qp];
  const Real a2 = (*_a2_prop)[_qp];
  const Real a3 = (*_a3_prop)[_qp];
  const Real p = (*_p_prop)[_qp];
  const Real d = _d[_qp];
  const Real eta = _eta;

  // degradation function
  _g[_qp] = std::pow((1-d),p)/(std::pow(1-d,p)+a1*d*(1+a2*d+a2*a3*std::pow(d,2)))*(1-_eta)+_eta;

  // Derivatives using quotient rule
  // U = (1-d)^p, V = a1*(d + a2*d^2 + a2*a3*d^3), D = U + V
  Real U   = std::pow(1-d, p);
  Real Up  = -p * std::pow(1-d, p-1);
  Real Up2 =  p*(p-1) * std::pow(1-d, p-2);

  Real V   = a1 * (d + a2*d*d + a2*a3*d*d*d);
  Real Vp  = a1 * (1 + 2*a2*d     + 3*a2*a3*d*d);
  Real Vpp = a1 * (    2*a2       + 6*a2*a3*d     );

  Real D   = U + V;
  Real Dp  = Up + Vp;
  Real Dpp = Up2 + Vpp;

  // first derivative g'
  Real N1  = Up*D - U*Dp;
  Real g0p = N1/(D*D);
  Real dg  = g0p * (1-eta);

  // second derivative g''
  Real N2  = Up2*D  - U*Dpp;
  Real g0pp = (N2*D - 2*N1*Dp)/(D*D*D);
  Real d2g  = g0pp * (1-eta);

  // store
  _dg_dd[_qp]    = dg;
  _d2g_dd2[_qp]  = d2g;
}
```

---

## 7. Phase-Field Evolution Equation and d=1 Behavior

### 7.1 Two Formulation Approaches

There are two main approaches to derive the phase-field evolution equation:

#### 7.1.1 Variational/Energy-Based Formulation (Current Implementation)

The current implementation uses energy minimization. The total energy functional:
```
Psi_total = integral [ g(d)*psi_e + Gc/c0*(alpha(d)/l + l*|grad(d)|^2) ] dV
```

The phase-field equation comes from the variational derivative (delta_Psi/delta_d = 0):
```
g'(d)*psie_active + (Gc/c0/l)*alpha'(d) - (2*Gc*l/c0)*Laplacian(d) = 0
```

In MOOSE, this is implemented as:
```
[Kernels]
  [diff]
    type = ADPFFDiffusion      # Handles: -2*Gc*l/c0 * Laplacian(d)
  []
  [source]
    type = ADPFFSource         # Handles: d(psi)/d(d) where psi = alpha*Gc/c0/l + g*psie_active
    free_energy = psi
  []
[]
```

#### 7.1.2 Driving Force Formulation (Miehe's Approach)

Miehe et al. use a different structure with explicit driving force (Eq. 13 in CMAME 2016):
```
eta * d_dot = (1-d)*H - [d - l^2*Laplacian(d)]
              --------   ---------------------
              driving    geometric resistance
              force
```

Key feature: The `(1-d)` factor is **explicitly** in front of the driving force H.

### 7.2 Behavior at d=1 (Fully Damaged Region)

A critical question: Does the phase-field equation properly stop evolving when d=1?

#### 7.2.1 AT1 Model Analysis

For AT1:
```
alpha(d) = d           -->  alpha'(d) = 1           (does NOT vanish at d=1!)
g(d) = (1-d)^2         -->  g'(d) = -2(1-d)         (vanishes at d=1)
```

At d=1, the equilibrium equation becomes:
```
0 + Gc/c0/l * 1 = 2*Gc*l/c0 * Laplacian(d)
     ---------
     Non-zero!
```

**Problem**: The `alpha'(d) = 1` term does not vanish, which could cause numerical issues.

#### 7.2.2 PF-CZM Model Analysis

For PF-CZM:
```
alpha(d) = 2d - d^2    -->  alpha'(d) = 2 - 2d = 2(1-d)    (VANISHES at d=1!)
g(d) = rational form   -->  g'(d) contains (1-d)^(p-1)     (vanishes at d=1 for p>=1)
```

At d=1, the equilibrium equation becomes:
```
0 + Gc/c0/l * 0 = 2*Gc*l/c0 * Laplacian(d)
```

**Both driving terms vanish!** This is a key advantage of PF-CZM.

#### 7.2.3 Miehe's Formulation

In Miehe's formulation:
```
eta * d_dot = (1-d)*H - [d - l^2*Laplacian(d)]
```

At d=1:
```
(1-d)*H = 0 * H = 0    (driving force explicitly vanishes!)
```

**Explicitly designed** to have zero driving force at d=1.

### 7.3 Comparison Summary: d=1 Behavior

| Model | alpha'(1) | g'(1) | Driving force at d=1 | Properly stops? |
|-------|-----------|-------|---------------------|-----------------|
| AT1 | 1 (non-zero!) | 0 | Gc/c0/l (non-zero) | **Problematic** |
| PF-CZM | 0 | 0 | 0 | **Yes** |
| Miehe | N/A | N/A | (1-d)*H = 0 | **Yes (explicit)** |

### 7.4 Physical Interpretation

When d = 1 (fully damaged):
1. **Stress state**: sigma = g(1)*sigma_eff ≈ 0 (traction-free)
2. **Driving force**: Should be zero to prevent further evolution
3. **Behavior**: Acts as a traction-free crack surface

For repeated/pulse loading:
```
                    Pulse Loading
                         |
                         v
    +--------------------+--------------------+
    |                    |                    |
    |   d ≈ 0           |||  d = 1  |||      d ≈ 0    |
    |   (intact)        ||| (crack) |||     (intact)  |
    |                   |||  band   |||               |
    |   carries         |||   no    |||     carries   |
    |   stress          ||| stress  |||     stress    |
    |                   |||         |||               |
    |                   |+=========+|                |
    |                    ^  crack  ^                  |
    |                    |  front  |                  |
    |              (damage grows here)                |
    +------------------------------------------------+
```

- **d=1 region**: Cannot accumulate more damage, acts as open crack
- **Crack front (d<1)**: Stress concentration, damage evolves, crack propagates
- **Intact region (d≈0)**: Carries load, may initiate new damage if threshold exceeded

### 7.5 Why PF-CZM Choice of alpha(d) = 2d - d^2 is Optimal

The choice `alpha(d) = 2d - d^2` gives `alpha'(d) = 2(1-d)` which:

1. **Vanishes at d=1**: Ensures no spurious driving force in fully damaged regions
2. **Equals 2 at d=0**: Provides proper initial damage evolution rate (xi = 2)
3. **Gives finite support**: Half bandwidth Du = pi*l/2 (bounded damage zone)
4. **Normalization c0 = pi**: Clean mathematical properties

This is why Wu (JMPS 2017) calls it the "optimal" crack geometric function for quasi-brittle materials.

---

## 8. References

1. Wu, J.Y. (2017). "A unified phase-field theory for the mechanics of damage and quasi-brittle failure." Journal of the Mechanics and Physics of Solids, 103, 72-99.

2. Wu, J.Y., & Nguyen, V.P. (2018). "A length scale insensitive phase-field damage model for brittle fracture." Journal of the Mechanics and Physics of Solids, 119, 20-42.

3. Miehe, C., & Mauthe, S. (2016). "Phase field modeling of fracture in multi-physics problems. Part III." Computer Methods in Applied Mechanics and Engineering, 304, 619-655.

4. RACCOON Documentation - Phase-Field Fracture Tutorials

5. Gupta et al. (2022). "An adaptive mesh refinement algorithm for phase-field fracture models: Application to brittle, cohesive, and dynamic fracture."

---

## 9. File Structure

```
code_pfczm/
├── PFCZM_Implementation_Plan.md     (this document)
├── elasticityhf_pfczm.i             (main elasticity input - dynamic solve)
├── fracturehf_pfczm.i               (fracture sub-app input - dynamic solve)
├── elasticityhf_static_pfczm.i      (main elasticity input - static solve)
└── fracturehf_static_pfczm.i        (fracture sub-app input - static solve)
```

### 9.1 Running the Simulation

**Static Solve (establish initial conditions):**
```bash
cd code_pfczm
mpirun -np 4 farms-opt -i elasticityhf_static_pfczm.i
```

**Dynamic Solve (with pulse loading):**
```bash
cd code_pfczm
mpirun -np 4 farms-opt -i elasticityhf_pfczm.i
```

---

## Appendix A: Wu's Paper Key Equations Reference

### A.1 Crack Surface Density (Eq. 2.5)
```
γ(d,∇d) = (1/c0)[α(d)/l + l|∇d|²]
c0 = 4∫₀¹ √α(β)dβ
```

### A.2 Stored Energy (Eq. 2.17)
```
ψ(ε,d) = ω(d)·ψ₀(ε)
ω(d) = (1-d)^p / [(1-d)^p + Q(d)]
```

### A.3 Failure Strength (Eq. 3.5)
```
ft = √(2·E·Gc·ξ / (c0·l·a1))
```

### A.4 Initial Slope (Eq. 3.7)
```
k0 = -(c0/4π)·(ft²/Gc)·[ξ(a2+p+1)-1]^(3/2) / ξ²
```

### A.5 Ultimate Crack Opening (Eq. 3.8)
```
For p = 2: wc = (2π·Gc)/(c0·ft)·√(ξ·P(1))
where P(1) = 1 + a2 + a2·a3
```

---

## Appendix B: Parameter Values for Current Implementation

### B.1 Material Properties
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| E | 30e9 | Pa | Young's modulus |
| nu | 0.3 | - | Poisson's ratio |
| Gc_const | 100 | N/m | Fracture energy |
| ft | 3e6 | Pa | Tensile strength |
| l | 2e-2 | m | Length scale |

### B.2 Derived PF-CZM Parameters
| Parameter | Expression | Value | Description |
|-----------|------------|-------|-------------|
| lch | E*Gc/ft² | 0.333 m | Irwin's characteristic length |
| c0 | π | 3.14159 | Normalization constant |
| xi | 2 | 2 | Initial slope of α(d) |
| a1 | (4/π)*lch/l | 21.22 | Degradation parameter |
| a2 | -0.5 | -0.5 | Linear softening |
| a3 | 0 | 0 | Linear softening |
| p | 2 | 2 | Power exponent |
| psic | 3*Gc/(8*l) | 187.5 J/m³ | Critical fracture energy |
| eta | 1e-6 | - | Residual stiffness |
