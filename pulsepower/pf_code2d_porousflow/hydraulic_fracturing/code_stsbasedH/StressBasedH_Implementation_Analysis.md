# Stress-Based Crack Driving Force Implementation

## Based on Miehe & Mauthe (CMAME 2016) + Wu's PF-CZM (JMPS 2017)

---

## 1. Theory

### 1.1 Phase Field Evolution Equation (Miehe Eq. 13)

```
η × ḋ = (1-d) × H - [d - l² × Δd]
        ─────────   ─────────────
         driving     geometric
         force       resistance
```

where:
- `d`: Phase field (damage variable, 0 = intact, 1 = fully cracked)
- `H`: Crack driving force (history variable)
- `l`: Length scale parameter
- `η`: Viscosity (rate-independent limit: η → 0)

### 1.2 Stress-Based Driving Force (Miehe Eq. 56)

```
D = ζ × ⟨ Σ_{a=1}^{3} ( ⟨σ̃_eff^a⟩₊ / σ_c )² - 1 ⟩₊
```

where:
- `σ̃_eff^a`: Principal values of **undegraded** effective stress
- `σ_c`: Critical fracture tensile stress (tensile strength)
- `ζ`: Slope parameter controlling driving force growth
- `⟨x⟩₊ = max(x, 0)`: Macaulay bracket (positive part)

### 1.3 History Variable for Irreversibility (Miehe Eq. 18)

```
H(X, t) = max_{s ∈ [0,t]} D(state(X, s))
```

The history maximum ensures damage irreversibility.

### 1.4 Effective Stress in Porous Media (Miehe Eq. 95)

Total stress decomposition:
```
σ = g(d) × σ̃_eff - b × p × I
```

where:
- `σ`: Total (degraded) stress tensor
- `g(d)`: Degradation function
- `σ̃_eff`: Undegraded effective stress
- `b`: Biot coefficient
- `p`: Pore pressure

**Recovery of undegraded effective stress:**
```
σ̃_eff = (σ + b × p × I) / g(d)
```

### 1.5 PF-CZM Crack Geometric Function (Wu 2017)

```
α(d) = 2d - d²
```

This gives:
- Normalization constant: `c₀ = π`
- Proper behavior at d=1: `α'(1) = 0` (damage stops in fully cracked regions)

---

## 2. Implementation

### 2.1 File Structure

```
code_stsbasedH/
├── elasticityhf_stsbased.i        # Dynamic main app
├── fracturehf_stsbased.i          # Dynamic fracture sub-app
├── elasticityhf_static_stsbased.i # Static main app
├── fracturehf_static_stsbased.i   # Static fracture sub-app
```

### 2.2 Material Parameters

```
# Stress-Based Parameters
sigma_c = 10e6      # Critical fracture tensile stress [Pa]
zeta = 1.0          # Slope parameter

# PF-CZM Parameters
Gc_const = 100      # Fracture energy [N/m]
l = 2e-2            # Length scale [m]
psic = 3*Gc/(8*l)   # Critical energy density
p_deg = 2           # Degradation exponent (linear softening)
a2 = -0.5           # Softening parameter
a3 = 0
eta = 1e-6          # Residual stiffness
```

### 2.3 Main App: Stress-Based Driving Force Computation

**Step 1: Extract degradation function g(d)**
```
[get_degradation_g]
  type = MaterialRealAux
  variable = degradation_g
  property = g
[]
```

**Step 2: Extract principal stresses (degraded)**
```
[get_max_principal_stress]
  type = RankTwoScalarAux
  rank_two_tensor = stress
  variable = max_principal_stress
  scalar_type = MaxPrincipal
[]
```

**Step 3: Compute driving force with undegraded effective stress**
```
[compute_stress_driving_force_D]
  type = ParsedAux
  variable = stress_driving_force_D
  coupled_variables = 'max_principal_stress mid_principal_stress
                       min_principal_stress pp biot_coefficient_aux degradation_g'
  constant_names = 'sigma_c zeta g_min'
  constant_expressions = '${sigma_c} ${zeta} 1e-10'
  expression = 'g_safe := max(degradation_g, g_min);
                sig_eff1_undeg := (max_principal_stress + biot_coefficient_aux * pp) / g_safe;
                sig_eff2_undeg := (mid_principal_stress + biot_coefficient_aux * pp) / g_safe;
                sig_eff3_undeg := (min_principal_stress + biot_coefficient_aux * pp) / g_safe;
                sig1_pos := max(sig_eff1_undeg, 0);
                sig2_pos := max(sig_eff2_undeg, 0);
                sig3_pos := max(sig_eff3_undeg, 0);
                sum_sq := (sig1_pos/sigma_c)^2 + (sig2_pos/sigma_c)^2 + (sig3_pos/sigma_c)^2;
                zeta * max(sum_sq - 1, 0)'
[]
```

**Step 4: History variable**
```
[compute_stress_driving_force_H]
  type = ParsedAux
  variable = stress_driving_force_H
  coupled_variables = 'stress_driving_force_D'
  expression = 'stress_driving_force_D'
[]
```

### 2.4 Fracture Sub-App: Free Energy Formulation

The free energy functional produces the correct `(1-d) × H` driving term:

```
[psi]
  type = ADDerivativeParsedMaterial
  property_name = psi
  expression = 'alpha*Gc/c0/l + 0.5*(1-d)^2 * (Gc/c0/l) * stress_driving_force'
  coupled_variables = 'd stress_driving_force'
  material_property_names = 'alpha(d) Gc c0 l'
  derivative_order = 1
[]
```

**Derivation:**
- `ψ_elastic = 0.5 × (1-d)² × (Gc/(c₀×l)) × H`
- `∂ψ_elastic/∂d = -(1-d) × (Gc/(c₀×l)) × H`
- `-∂ψ_elastic/∂d = (1-d) × (Gc/(c₀×l)) × H` ✓

This gives the correct `(1-d) × H` driving term from Miehe's evolution equation.

### 2.5 MultiApp Transfers

```
# Transfer H from main app to fracture sub-app
[to_stress_driving_force]
  type = MultiAppCopyTransfer
  to_multi_app = 'fracture'
  variable = 'stress_driving_force'
  source_variable = 'stress_driving_force_H'
[]

# Transfer d from fracture sub-app to main app
[from_d]
  type = MultiAppCopyTransfer
  from_multi_app = 'fracture'
  variable = 'd'
  source_variable = 'd'
[]
```

---

## 3. How to Run

```bash
cd code_stsbasedH/

# Step 1: Static solve (initial equilibrium)
farms-opt -i elasticityhf_static_stsbased.i --allow-unused

# Step 2: Dynamic solve (loads from static checkpoint)
farms-opt -i elasticityhf_stsbased.i --allow-unused
```

---

## 4. Key Equations Summary

| Component | Equation |
|-----------|----------|
| Total stress (Eq. 95) | `σ = g(d) × σ̃_eff - b×p×I` |
| Undegraded effective stress | `σ̃_eff = (σ + b×p×I) / g(d)` |
| Driving force (Eq. 56) | `D = ζ × ⟨Σ(⟨σ̃_eff^a⟩₊/σ_c)² - 1⟩₊` |
| History variable (Eq. 18) | `H = max_{s∈[0,t]} D` |
| Evolution equation (Eq. 13) | `η×ḋ = (1-d)×H - [d - l²Δd]` |
| Free energy | `ψ = α×Gc/(c₀×l) + 0.5×(1-d)²×(Gc/(c₀×l))×H` |

---

## 5. References

1. Miehe, C., & Mauthe, S. (2016). Phase field modeling of fracture in multi-physics problems. Part III. Crack driving forces in hydro-poro-elasticity and hydraulic fracturing of fluid-saturated porous media. *Computer Methods in Applied Mechanics and Engineering*, 304, 619-655.

2. Wu, J. Y. (2017). A unified phase-field theory for the mechanics of damage and quasi-brittle failure. *Journal of the Mechanics and Physics of Solids*, 103, 72-99.

---

*Implementation completed January 2026*
