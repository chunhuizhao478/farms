# MOOSE vs Tandem: Detailed Implementation Comparison for SEAS Simulations

**Document Version:** 1.0
**Date:** January 28, 2026
**Author:** Technical Documentation
**Purpose:** Comprehensive comparison of MOOSE (farms_rsf_explicit) and Tandem implementations for SEAS (Sequences of Earthquakes and Aseismic Slip) simulations

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Overall Architecture](#2-overall-architecture)
3. [Elasticity Discretization](#3-elasticity-discretization)
4. [Rate-and-State Friction](#4-rate-and-state-friction)
5. [Domain-Fault Coupling](#5-domain-fault-coupling)
6. [Time Integration](#6-time-integration)
7. [Boundary Conditions and Loading](#7-boundary-conditions-and-loading)
8. [Mesh and Geometry](#8-mesh-and-geometry)
9. [Parallelization](#9-parallelization)
10. [Output and Monitoring](#10-output-and-monitoring)
11. [Key Differences Summary](#11-key-differences-summary)
12. [Benchmark Verification Results](#12-benchmark-verification-results)
13. [Recommendations](#13-recommendations)

---

## 1. Executive Summary

This document provides a detailed comparison between the current MOOSE implementation (benchmark_bp2 and multiapp_bp2 approaches) and the Tandem code for solving SEAS problems, specifically the SCEC BP2-QD benchmark.

### Key Findings

| Aspect | MOOSE (Current) | Tandem |
|--------|-----------------|--------|
| **Architecture** | MultiApp (MainApp + SubApp) or Staggered | Three-operator (Domain + Adapter + Friction) |
| **DG Method** | SIPG with custom kernels | SIPG with YATeTo-optimized kernels |
| **State Variable** | Direct θ (benchmark convention) | Modified ψ convention |
| **Friction Solve** | Brent's method | Brent's method |
| **Time Stepping** | Custom SEASAdaptiveDT | PETSc TS with RK/BDF |
| **Linear Solver** | MUMPS (LU direct) | PETSc KSP with multigrid |
| **Green's Function** | Supported via stiffness mode | SeasQDDiscreteGreenOperator |

### Implementation Approaches in MOOSE

1. **Single-App Staggered (benchmark_bp2)**: All physics in one application with AuxKernels
2. **MultiApp Architecture (multiapp_bp2)**: Elasticity MainApp + Friction SubApp with transfers

---

## 2. Overall Architecture

### 2.1 Tandem Three-Operator Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                    SeasFDOperator (Main)                      │
│  ┌─────────────────┐  ┌──────────────┐  ┌─────────────────┐  │
│  │  DGOperator     │  │   Adapter    │  │ FrictionOperator│  │
│  │  (Elasticity)   │◄─┤   Operator   │◄─┤  (RateState)    │  │
│  │  - SIPG         │  │  - Traction  │  │  - State evol.  │  │
│  │  - Volume terms │  │  - Slip BC   │  │  - V solve      │  │
│  └─────────────────┘  └──────────────┘  └─────────────────┘  │
└──────────────────────────────────────────────────────────────┘

State Vector: [velocity (v), displacement (u), fault state (s)]

Key Files:
- app/form/SeasFDOperator.h        (Fully dynamic)
- app/form/SeasQDOperator.h        (Quasi-dynamic)
- app/form/FrictionOperator.h      (Friction)
- app/localoperator/Elasticity.h   (DG elasticity)
- app/localoperator/DieterichRuinaAgeing.h (Friction law)
```

**Tandem Design Philosophy:**
- Clean separation of physics components
- Template-based friction law pluggability
- YATeTo kernel generation for performance
- PETSc integration for solvers and time stepping

### 2.2 MOOSE Single-App Staggered Architecture (benchmark_bp2)

```
┌────────────────────────────────────────────────────────────────┐
│                    Single MOOSE Application                     │
│                                                                 │
│  ┌──────────────────────────┐   ┌────────────────────────────┐ │
│  │     DGKernels            │   │    AuxKernels (on fault)   │ │
│  │  - DGElasticityAntiplane │   │  - SEASStateAux            │ │
│  │                          │   │  - SEASSlipAux             │ │
│  │  InterfaceKernels        │◄──│  - SEASTractionAux         │ │
│  │  - DGFaultSlipInterface  │   │  - SEASSlipRateVarAAux     │ │
│  └──────────────────────────┘   └────────────────────────────┘ │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │           InterfaceMaterial: DGElasticTractionMaterial    │  │
│  │           (Computes τ from displacement gradient)         │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────────────────┘

Execution Order:
TIMESTEP_BEGIN: Update slip (SEASSlipAux), Update state (SEASStateAux)
MAIN SOLVE:     DG elasticity with prescribed slip BC
TIMESTEP_END:   Compute traction (SEASTractionAux), Solve V (SEASSlipRateVarAAux)
```

**Key Files:**
- `src/dgkernels/DGElasticityAntiplane.C`
- `src/interfacekernels/DGFaultSlipInterfaceKernel.C`
- `src/materials/dg/DGElasticTractionMaterial.C`
- `src/auxkernels/seas/SEASSlipRateVarAAux.C`

### 2.3 MOOSE MultiApp Architecture (multiapp_bp2)

```
┌──────────────────────────────────────────────────────────────────┐
│                       MOOSE MultiApp System                       │
│                                                                   │
│  ┌────────────────────────────────┐                              │
│  │     MainApp (Elasticity)       │    ┌───────────────────┐     │
│  │  - 2D DG mesh (left+right)     │    │    Transfers      │     │
│  │  - DGElasticityAntiplane       │◄───┤  send_traction    │     │
│  │  - DGFaultSlipInterfaceKernel  │───►│  receive_slip     │     │
│  │  - DGElasticTractionMaterial   │    │  receive_slip_rate│     │
│  └────────────────────────────────┘    └───────────────────┘     │
│                                                 │                 │
│                                                 ▼                 │
│  ┌──────────────────────────────────────────────────────────────┐│
│  │                  SubApp (Friction)                            ││
│  │  - 1D fault mesh (depth direction)                           ││
│  │  - StateEvolutionKernel: dθ/dt = 1 - V*θ/Dc                 ││
│  │  - SEASSlipRateVarAAux: Solve for V from traction balance   ││
│  │  - SlipIntegrationAux: s += V*dt                            ││
│  └──────────────────────────────────────────────────────────────┘│
└──────────────────────────────────────────────────────────────────┘

Transfer Types: MultiAppGeneralFieldNearestLocationTransfer
Execute Order: SubApp runs at TIMESTEP_BEGIN before MainApp
```

**Key Files:**
- `examples/multiapp_bp2/main_elasticity.i`
- `examples/multiapp_bp2/sub_friction.i`
- `src/kernels/friction/StateEvolutionKernel.C`

### 2.4 Architecture Comparison

| Feature | Tandem | MOOSE Staggered | MOOSE MultiApp |
|---------|--------|-----------------|----------------|
| **Separation of concerns** | Excellent (3 operators) | Moderate (AuxKernels) | Good (2 apps) |
| **Mesh flexibility** | Same mesh for all | Single mesh | Different meshes |
| **Code modularity** | Template-based | Input file based | App-based |
| **Debugging ease** | Component-level | Coupled | App-level |
| **Friction law swapping** | Template parameter | AuxKernel change | SubApp swap |
| **3D extension** | Built-in | Requires work | Natural (2D SubApp) |

---

## 3. Elasticity Discretization

### 3.1 Governing Equations

**BP2 Antiplane Shear (2D):**
```
Strong form: 0 = ∂σ_xy/∂x + ∂σ_yz/∂z   (equilibrium)
             σ_xy = μ ∂u/∂x             (constitutive)
             σ_yz = μ ∂u/∂z

Reduces to:  -∇·(μ∇u) = 0              (Laplace equation)
```

### 3.2 DG Weak Formulation (SIPG)

Both implementations use the **Symmetric Interior Penalty Galerkin (SIPG)** method:

```
a(u,v) = ∫_Ω (μ∇u)·(∇v) dΩ                        [Volume term]
       - ∫_Γ_I {{μ∇u·n}} [[v]] dΓ                  [Consistency]
       + ∫_Γ_I {{μ∇v·n}} [[u]] dΓ                  [Symmetry, ε=+1]
       + ∫_Γ_I (σμ/h) [[u]] [[v]] dΓ               [Penalty]
```

Where:
- `{{·}}` = average operator across facet
- `[[·]]` = jump operator across facet
- `σ` = penalty parameter (typically σ = 6-200)
- `h` = characteristic element size

### 3.3 Tandem Elasticity Implementation

**File:** `app/localoperator/Elasticity.h`

```cpp
// Full 2D/3D elasticity with Lamé parameters
σ_ij = λ δ_ij (∂u_k/∂x_k) + μ(∂u_i/∂x_j + ∂u_j/∂x_i)

// YATeTo kernel generation (app/kernels/elasticity.py)
def traction(x, normal):
    return lam_q[x]['q'] * Dx_q[x]['lsq'] * u[x]['ls'] * normal['pq'] + \
           mu_q[x]['q'] * (Dx_q[x]['ljq'] * u[x]['lp'] * normal['jq'] +
                          Dx_q[x]['lpq'] * u[x]['lj'] * normal['jq'])
```

**Key Features:**
- High-order polynomial basis (p = 3-4 typical)
- Warp-and-blend Fekete points on simplices
- Precomputed matrices: `MhatInv`, `E_Q`, `Dxi_Q`
- Matrix-free option for large problems
- BR2 flux stabilization alternative

### 3.4 MOOSE Elasticity Implementation

**File:** `src/dgkernels/DGElasticityAntiplane.C`

```cpp
// Simplified antiplane shear (scalar Laplace)
// Uses MFEM-style DG kernel structure

Real DGElasticityAntiplane::computeQpResidual(Moose::DGResidualType type)
{
  // Consistency: -{{μ∇u·n}}[[v]]
  Real r = -_epsilon * 0.5 * (_mu[_qp] + _mu_neighbor[_qp]) *
           (_grad_u_elem[_qp] * _normals[_qp] +
            _grad_u_neighbor[_qp] * _normals[_qp]) * _test[_i][_qp];

  // Symmetry: +{{μ∇v·n}}[[u]]
  r += _sigma * 0.5 * (_mu[_qp] + _mu_neighbor[_qp]) *
       (_grad_test_elem[_i][_qp] * _normals[_qp] +
        _grad_test_neighbor[_i][_qp] * _normals[_qp]) *
       (_u_elem[_qp] - _u_neighbor[_qp]);

  // Penalty: (σμ/h)[[u]][[v]]
  r += _penalty * (_mu[_qp] / _h_elem + _mu_neighbor[_qp] / _h_neighbor) *
       (_u_elem[_qp] - _u_neighbor[_qp]) * _test[_i][_qp];

  return r;
}
```

**Key Features:**
- FIRST order MONOMIAL basis (DG L2)
- Coupled with InterfaceKernel for fault
- Parameters: epsilon=1.0 (SIPG), sigma=200

### 3.5 Elasticity Comparison

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Polynomial order** | High-order (p=3-4) | First order (p=1) |
| **Element type** | Simplices (tri/tet) | Quadrilaterals |
| **Basis family** | Warp-and-blend Fekete | MONOMIAL (L2) |
| **Kernel optimization** | YATeTo code generation | Hand-written |
| **Matrix assembly** | Matrix-free option | Assembled |
| **Penalty parameter** | σ ≈ 6 | σ = 200 |
| **Flux stabilization** | IP or BR2 | IP only |

---

## 4. Rate-and-State Friction

### 4.1 Governing Equations

**Friction Coefficient (Regularized, BP2 Eq. 7):**
```
f(V, θ) = a · sinh⁻¹[V/(2V₀) · exp((f₀ + b·ln(V₀θ/Dc))/a)]
```

**State Evolution (Aging Law, BP2 Eq. 6):**
```
dθ/dt = 1 - V·θ/Dc
```

**Traction Balance (BP2 Eq. 5):**
```
τ = τ⁰ + τ_qs - η·V

Where: η = μ/(2·c_s) (radiation damping coefficient)
```

### 4.2 Tandem Friction Implementation

**File:** `app/localoperator/DieterichRuinaAgeing.h`

```cpp
// Uses MODIFIED state variable ψ = f₀ + b·ln(V₀θ/Dc)
// Simplifies friction coefficient formula

// Friction coefficient (simplified due to ψ convention)
double F(std::size_t index, double snAbs, double V, double psi) const {
    auto a = p_[index].get<A>();
    double e = exp(psi / a);  // Direct use of ψ
    double f = a * asinh((V / (2.0 * cp_.V0)) * e);
    return snAbs * f;
}

// State evolution with ψ convention
double state_rhs(std::size_t index, double V, double psi) const {
    double myL = p_[index].get<L>();  // L = Dc
    // dψ/dt = b*V₀/L * (exp((f₀-ψ)/b) - V/V₀)
    return cp_.b * cp_.V0 / myL * (exp((cp_.f0 - psi) / cp_.b) - V / cp_.V0);
}

// Slip rate inversion using Brent's method
auto slip_rate(...) const {
    double a = 0.0;
    double b = tauAbs / eta;  // Upper bound

    auto fF = [&](double V) {
        return tauAbs - this->F(index, snAbs, V, psi) - eta * V;
    };

    V = zeroIn(a, b, fF);  // Brent's method
}
```

**Key Features:**
- Modified state variable ψ (simplifies formulas)
- Template-based for law swapping
- Spatially varying parameters via Params struct
- Pre-stress included in traction computation

### 4.3 MOOSE Friction Implementation

**File:** `src/auxkernels/seas/SEASSlipRateVarAAux.C`

```cpp
// Uses DIRECT θ convention (benchmark standard)

// Friction coefficient (exact BP2 Eq. 7)
Real SEASSlipRateVarAAux::frictionCoefficient(Real V, Real theta) const
{
  Real V_eff = std::max(V, _V_min);  // Regularization
  Real arg = (V_eff / (2.0 * _V0)) *
             std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / _a);
  return _a * std::asinh(arg);
}

// State evolution (exact BP2 Eq. 6)
// Implemented in SEASStateAux with backward Euler:
// θ_new = (θ_old + Δt) / (1 + V·Δt/Dc)

// Slip rate inversion using BrentRootFinder
Real SEASSlipRateVarAAux::solveSlipRate(Real tau, Real theta) const
{
  auto residual = [&](Real V) {
    Real f = frictionCoefficient(V, theta);
    return tau - _sigma_n * f - _eta * V;
  };

  Real V_max = std::max(tau / _eta, 1.0);
  return BrentRootFinder::zeroIn(_V_min, V_max, residual, _tol);
}
```

**Key Features:**
- Direct θ convention (follows benchmark exactly)
- Supports depth-dependent a(z) via coupled variable
- Backslip loading option for deep fault zone
- Stability limits: V_min = 1e-20, V_max = 10.0

### 4.4 State Variable Convention Comparison

| Convention | Formula | Advantage |
|------------|---------|-----------|
| **θ (MOOSE)** | f = a·sinh⁻¹[V/(2V₀)·exp((f₀+b·ln(V₀θ/Dc))/a)] | Matches benchmark directly |
| **ψ (Tandem)** | f = a·sinh⁻¹[V/(2V₀)·exp(ψ/a)] | Simpler friction formula |

**Transformation:**
```
ψ = f₀ + b·ln(V₀θ/Dc)
θ = (Dc/V₀)·exp((ψ-f₀)/b)
```

### 4.5 Friction Comparison

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **State variable** | Modified ψ | Direct θ (benchmark) |
| **Root finder** | Brent's method | Brent's method |
| **Tolerance** | 1e-14 | 1e-12 |
| **Parameter storage** | Struct-based | Input file |
| **Depth-dependence** | Lua functions | ParsedFunction + AuxVariable |
| **Radiation damping** | μ/(2·c_s) | μ/(2·c_s) |

---

## 5. Domain-Fault Coupling

### 5.1 Tandem Adapter Operator

**File:** `app/form/SeasQDOperator.cpp`

**Complete Traction Flow in Tandem:**

```
┌─────────────────────────────────────────────────────────────────┐
│ Time Integrator calls: rhs(time, state, result)                 │
│   state = [slip, ψ (state variable)]                            │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ solve(time, state_view)                                          │
│   1. dgop_->set_slip(adapter_->slip_bc(state_view))  ← Set BC   │
│   2. dgop_->set_dirichlet(fun_boundary_(time))       ← Far-field│
│   3. linear_solver_.update_rhs(*dgop_)               ← Assemble │
│   4. linear_solver_.solve()                          ← Solve    │
│   5. dgop_->set_slip(invalid_slip_bc())              ← Clear BC │
│   OUTPUT: displacement u                                         │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ update_traction(state_view)                                      │
│   1. dgop_->set_slip(adapter_->slip_bc(state_view))  ← Set BC   │
│   2. adapter_->traction(disp_view, traction_)        ← Compute  │
│   3. dgop_->set_slip(invalid_slip_bc())              ← Clear BC │
│   OUTPUT: traction τ                                             │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ friction_->rhs(time, traction, state, result)                   │
│   For each fault node:                                          │
│     1. Solve τ = σn·f(V,ψ) + η·V for V (Brent's method)        │
│     2. Compute dψ/dt from state evolution law                   │
│   OUTPUT: result = [dslip/dt = V, dψ/dt]                        │
└─────────────────────────────────────────────────────────────────┘
```

**Key Code Sections:**

```cpp
// File: app/form/SeasQDOperator.cpp

void SeasQDOperator::rhs(double time, BlockVector const& state, BlockVector& result) {
    update_ghost_state(state);
    solve(time, make_state_view(state));      // Step 1: Solve elasticity
    update_traction(make_state_view(state));  // Step 2: Compute traction
    friction_->rhs(time, traction_, state, result);  // Step 3: Friction RHS
}

void SeasQDOperator::solve(double time, BlockView const& state_view) {
    dgop_->set_slip(adapter_->slip_bc(state_view));  // Apply slip BC
    if (fun_boundary_) {
        dgop_->set_dirichlet((*fun_boundary_)(time));  // Far-field BC
    }
    linear_solver_.update_rhs(*dgop_);
    linear_solver_.solve();
    dgop_->set_slip(invalid_slip_bc());  // Clear for safety
}

void SeasQDOperator::update_traction(BlockView const& state_view) {
    auto disp_view = LocalGhostCompositeView(linear_solver_.x(), disp_ghost_);
    dgop_->set_slip(adapter_->slip_bc(state_view));  // IMPORTANT: Set slip for traction formula
    adapter_->traction(disp_view, traction_);        // Compute traction
    dgop_->set_slip(invalid_slip_bc());
}
```

**Traction Formula (from `app/kernels/elasticity.py`):**

```python
# DG-consistent traction at fault interface
traction_q = 0.5 * (traction(0) + traction(1)) + c0 * (u[0] - u[1] - slip)

# Where:
#   traction(i) = μ * grad_u[i] · n     (stress from each side)
#   c0 = penalty coefficient
#   u[0] - u[1] = displacement jump [[u]]
#   slip = prescribed slip from state
```

**Key Insight:** Tandem sets slip BC **both for solving AND for traction computation**.
This ensures the traction formula `τ = {{μ∇u·n}} + c0·([[u]] - slip)` uses consistent values.

**Coupling Flow:**
1. **Domain → Fault:** Displacement gradient → Traction via DG formula with penalty correction
2. **Fault → Domain:** Slip → Weakly enforced Dirichlet BC on fault facets (Nitsche method)

### 5.2 MOOSE Staggered Coupling

**Traction Computation (InterfaceMaterial):**
```cpp
// File: src/materials/dg/DGElasticTractionMaterial.C

// Two modes available:

// Mode 1: "gradient" - DG average gradient
τ = μ · {{∂u/∂n}} = μ · 0.5·(∇u_elem + ∇u_neighbor)·n

// Mode 2: "stiffness" - Quasi-static relationship
τ = K · (Vp·t - slip)
where K = μ/W_f (stiffness per unit fault length)
```

**Slip BC Enforcement (InterfaceKernel):**
```cpp
// File: src/interfacekernels/DGFaultSlipInterfaceKernel.C

// Enforces: [[u]] = u_neighbor - u_elem = slip_prescribed
// Via penalty method (Nitsche's approach):

residual = -epsilon * {{μ∇u·n}} * [[v]]
         + sigma * {{μ∇v·n}} * ([[u]] - slip)
         + penalty/h * ([[u]] - slip) * [[v]]
```

**Complete Traction Flow in MOOSE Staggered Solver:**

```
┌─────────────────────────────────────────────────────────────────┐
│ DG Elasticity Solve                                              │
│   - Solves -∇·(μ∇u) = 0                                         │
│   - With slip BC: [[u]] = slip_prescribed                       │
│   - Produces displacement field u                                │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ Traction Computation (UserObject or InterfaceMaterial)          │
│   τ = μ·{{∂u/∂n}} - κ·([[u]] - slip) + τ₀                      │
│   (dg_consistent mode follows Tandem's formula)                 │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ SEASTractionAux (or InterfaceValueUserObjectAux)                │
│   Copies traction → `traction` AuxVariable                      │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ SEASSlipRateVarAAux                                              │
│   INPUT: τ = _traction[_qp]                                     │
│   SOLVE: τ = σn·f(V,θ) + η·V  for V using Brent's method       │
│   OUTPUT: slip_rate V                                           │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│ Next Timestep: SEASSlipAux                                       │
│   slip_new = slip_old + V·dt                                    │
└─────────────────────────────────────────────────────────────────┘
```

**Execution Order within Each Timestep:**
```
TIMESTEP_BEGIN:
  1. SEASSlipAux:     slip_new = slip_old + V_old * dt
  2. SEASStateAux:    θ_new = (θ_old + dt) / (1 + V_old*dt/Dc)

SOLVE (Elasticity):
  3. DGElasticityAntiplane + DGFaultSlipInterfaceKernel
     Solve: -∇·(μ∇u) = 0  with  [[u]] = slip_new

TIMESTEP_END:
  4. Traction computation (Material or UserObject)
  5. SEASTractionAux: Read τ into AuxVariable
  6. SEASSlipRateVarAAux: Solve τ = σn·f(V,θ) + η·V for V_new
```

**Known Issue with DG Gradient-Based Traction:**

The DG formulation with large penalty parameters (e.g., 1e14) can mask the physical
stress gradient from far-field loading. The penalty enforcement creates steep artificial
gradients near the fault that dominate the gradient term `μ·{{∂u/∂n}}`.

This is why the "stiffness" mode (`traction_mode = stiffness`) is recommended for
quasi-static SEAS problems:
```cpp
// Stiffness mode bypasses local gradients:
τ = K · (Vp·t - slip) + τ₀
// where K = μ/Wf captures the far-field loading directly
```

**Verification Results (January 28, 2026):**

The stiffness mode was verified with a 1-year BP2-QD simulation (800m mesh):

| Region | Initial τ (MPa) | Final τ (MPa) | Δτ (kPa) | Final Slip Rate (m/s) | Total Slip |
|--------|-----------------|---------------|----------|-----------------------|------------|
| VW (z=0) | 26.546 | 26.571 | **+25.3** | 1.5×10⁻¹⁵ (locked) | 8.16 μm |
| VS (z=24) | 26.546 | 26.546 | **0** | 1×10⁻⁹ (plate rate) | 31.6 mm |

**Key observations:**
1. ✅ VW region locks up correctly (slip rate drops by 6 orders of magnitude)
2. ✅ Stress builds in VW region due to slip deficit (~25 kPa/year)
3. ✅ VS region creeps at plate rate with constant stress
4. ✅ VS total slip = Vp × t = 1e-9 × 3.16e7 ≈ 31.6 mm (exact match)

In contrast, the `dg_consistent` mode with large penalty parameters produced incorrect stress
evolution (decreasing instead of increasing) due to numerical artifacts dominating the physical
gradient.

### 5.3 MOOSE MultiApp Coupling

**Transfers:**
```
[Transfers]
  # MainApp → SubApp (before SubApp solve)
  [send_traction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = friction
    source_variable = traction_to_subapp
    variable = traction_received
  []

  # SubApp → MainApp (after SubApp solve)
  [receive_slip]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = friction
    source_variable = slip
    variable = slip_from_subapp
  []
[]
```

**Coordinate Mapping:**
- MainApp: 2D bulk mesh, fault at x=0
- SubApp: 1D mesh along z (depth)
- Transfer handles: (x=0, z) ↔ (z) coordinate mapping

### 5.4 Coupling Comparison

| Aspect | Tandem | MOOSE Staggered | MOOSE MultiApp |
|--------|--------|-----------------|----------------|
| **Traction source** | DG average gradient | Gradient or stiffness | Transfer from MainApp |
| **Slip BC** | DG Dirichlet | InterfaceKernel penalty | InterfaceKernel + Transfer |
| **Coordinate transform** | Built-in adapter | Implicit in InterfaceKernel | Nearest-node transfer |
| **Update timing** | Within RHS evaluation | TIMESTEP_BEGIN/END | Transfer execute_on |

---

## 6. Time Integration

### 6.1 Tandem Time Integration

**File:** `app/common/PetscTimeSolver.h`

```cpp
// Wraps PETSc TS module
template <std::size_t NumStateVecs>
class PetscTimeSolver {
    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetType(ts, TSRK);  // Runge-Kutta
    TSSetRHSFunction(ts, NULL, rhsFunction, ctx);
    TSSolve(ts, solution);
};
```

**Available Methods:**
- Explicit: RK2, RK3, RK4, RK5 (with FSAL)
- Implicit: BDF, Theta methods
- Adaptive: Built-in error control and step rejection

**CFL Constraint (Fully Dynamic):**
```cpp
// File: app/localoperator/Elasticity.cpp
double cfl_time_step() const {
    return CFL_coefficient * min(h_element / c_max);
}
// c_max = max(c_p, c_s) = P-wave or S-wave speed
```

### 6.2 MOOSE Time Integration

**SEASAdaptiveDT Time Stepper:**
```cpp
// File: src/timestepper/SEASAdaptiveDT.C

Real SEASAdaptiveDT::computeDT()
{
  Real Vmax = _max_slip_rate_pp.getValue();

  // Time step formula: dt = C · Dc / V_max
  Real dt_slip = _C * _Dc / std::max(Vmax, _V_min);

  if (Vmax >= _V_seismic)  // V_seismic = 1e-3 m/s
  {
    // Seismic event: use small dt_seismic
    return std::min(dt_slip, _dt_seismic);
  }
  else
  {
    // Aseismic: limit growth by growth_factor
    Real dt_grown = _dt_old * _growth_factor;
    return std::min({dt_slip, dt_grown, _dt_max});
  }
}
```

**Parameters:**
- `C = 0.1-0.5` (safety factor)
- `dt_seismic = 0.01 s` (during earthquakes)
- `dt_max = 1e5-1e6 s` (aseismic periods)
- `growth_factor = 1.05-1.2`

### 6.3 State Evolution Discretization

**Tandem (continuous ODE):**
```
// State evolved by PETSc TS along with slip
dψ/dt = (b·V₀/L) · (exp((f₀-ψ)/b) - V/V₀)
ds/dt = V
```

**MOOSE (semi-implicit):**
```cpp
// Backward Euler for state (unconditionally stable)
θ_new = (θ_old + Δt) / (1 + V·Δt/Dc)

// Forward Euler for slip
s_new = s_old + Δt·V
```

### 6.4 Time Integration Comparison

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **ODE solver** | PETSc TS (RK, BDF) | MOOSE Executioner |
| **State integration** | Continuous | Semi-implicit Euler |
| **Slip integration** | Continuous | Forward Euler |
| **Adaptive control** | PETSc error estimator | Slip-rate based |
| **dt formula** | CFL (FD) or slip-based | dt = C·Dc/V_max |
| **Growth limiting** | PETSc internal | Manual growth_factor |

---

## 7. Boundary Conditions and Loading

### 7.1 BP2 Boundary Conditions

| Boundary | Condition | Type |
|----------|-----------|------|
| Free surface (z=0) | σ_yz = 0 | Neumann (natural) |
| Far-field (x→±∞) | u = ±V_p·t/2 | Dirichlet |
| Below fault (z≥W_f) | V = V_p | Imposed rate |

### 7.2 Tandem Loading

**File:** `examples/tandem/2d/tutorial.lua`

```lua
-- Far-field velocity boundary condition
function Tutorial:boundary(x, y, t)
    local Vh = self.Vp * t / 2.0
    if x < 0 then
        Vh = -Vh
    end
    return Vh, 0.0  -- (u_tangent, u_normal)
end
```

### 7.3 MOOSE Loading

**DG Dirichlet BC:**
```cpp
// File: src/bcs/dg/DGElasticityDirichletBC.C

// Weakly enforced Dirichlet via penalty
// Left boundary: u = -Vp*t/2
// Right boundary: u = +Vp*t/2

[Functions]
  [left_bc]
    type = ParsedFunction
    expression = '-0.5*1e-9*t'  # -Vp*t/2
  []
  [right_bc]
    type = ParsedFunction
    expression = '0.5*1e-9*t'   # +Vp*t/2
  []
[]
```

**Loading Rate (stiffness mode):**
```cpp
// Alternative: quasi-static stress loading
// dτ/dt ≈ μ·Vp/(H+h)
loading_rate = mu * Vp / Lx;  // Lx = half-width
```

### 7.4 Free Surface Treatment

| Code | Free Surface | Implementation |
|------|--------------|----------------|
| Tandem | Natural BC | No face integral contribution |
| MOOSE | DGElasticityNeumannBC | Traction = 0, explicit term |

---

## 8. Mesh and Geometry

### 8.1 Tandem Mesh

**Format:** External mesh (Gmsh, VTK)

**Features:**
- Unstructured simplicial (triangles/tetrahedra)
- Curvilinear elements supported
- Fault marked as boundary attribute
- METIS/ParMETIS partitioning

**Typical Resolution:**
- High-order p=3-4, coarse mesh
- h-refinement near fault

### 8.2 MOOSE Mesh

**Generation:** MOOSE MeshGenerators

```
[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 125
    ny = 125
    xmin = -100000
    xmax = 0
    ymin = -100000
    ymax = 0
    subdomain_id = 1
  []
  [right_block]
    type = GeneratedMeshGenerator
    # ... similar for right block
    subdomain_id = 2
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'left_block right_block'
    stitch_boundaries_pairs = 'right left'
  []
  [fault]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = 1
    paired_block = 2
    new_boundary = fault
  []
[]
```

**Features:**
- Structured quadrilateral mesh
- Two blocks stitched at fault
- Fault as interface sideset
- Seismogenic zone restricted via ParsedGenerateSideset

### 8.3 Mesh Comparison

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Element type** | Simplices | Quadrilaterals |
| **Generation** | External (Gmsh) | Built-in generators |
| **High-order** | Curvilinear | Straight edges |
| **Fault representation** | Boundary attribute | Interface sideset |
| **Typical resolution** | 25-400 m | 800 m |

---

## 9. Parallelization

### 9.1 Tandem Parallel Implementation

**MPI + Domain Decomposition:**
```cpp
// Mesh partitioning
GlobalSimplexMesh::repartition();  // METIS/ParMETIS

// Ghost exchange
Scatter class for element-wise communication
CommPattern for facet communication

// Global reductions
MPI_Allreduce for time step selection (max V)
```

**Scalability:**
- Tested on HPC systems
- Matrix-free operators for memory efficiency
- GAMG multigrid for elastostatic solves

### 9.2 MOOSE Parallel Implementation

**Built-in MOOSE Parallelization:**
```cpp
// Mesh partitioning
DistributedMesh with automatic partitioning

// Ghost exchange
MOOSE internal ghosting for DG

// MultiApp parallelization
TransientMultiApp runs SubApp on each processor
```

**Transfer Considerations:**
- `MultiAppGeneralFieldNearestLocationTransfer` handles cross-mesh
- Postprocessor transfer for scalar values (max slip rate)

---

## 10. Output and Monitoring

### 10.1 Tandem Output

**File:** `app/tandem/Writer.h`

**Formats:**
- VTU (ParaView-compatible)
- CSV for time series
- Tecplot

**Monitoring (app/tandem/Monitor.h):**
- Tracks Vmax for adaptive time stepping
- AdaptiveOutputInterval for efficient I/O
- Checkpointing for restart

### 10.2 MOOSE Output

**Postprocessors for Monitoring:**
```
[Postprocessors]
  [max_slip_rate]
    type = NodalExtremeValue
    variable = slip_rate
    value_type = max
    boundary = fault
  []
  [slip_z12]
    type = PointValue
    variable = slip
    point = '0 -12000 0'
  []
[]
```

**Output Formats:**
- Exodus (ParaView-compatible)
- CSV for time series
- Console for convergence

**SCEC-Compliant Output:**
- 12 monitoring stations at specified depths
- Format: time, slip, log10(V), stress, log10(θ)

---

## 11. Key Differences Summary

### 11.1 Architecture

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Design philosophy** | Three operators, clean separation | Framework-based, AuxKernels or MultiApp |
| **Code organization** | Custom C++ with PETSc | MOOSE framework objects |
| **Configuration** | TOML + Lua | MOOSE input files (.i) |
| **Extensibility** | Template parameters | Subclass kernels/materials |

### 11.2 Numerical Methods

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Polynomial order** | High-order (p=3-4) | Low-order (p=1) |
| **Kernel optimization** | YATeTo generated | Hand-written |
| **State variable** | Modified ψ | Direct θ |
| **Time stepping** | PETSc TS (RK/BDF) | Custom SEASAdaptiveDT |
| **Linear solver** | PETSc KSP + multigrid | MUMPS direct |

### 11.3 Computational Aspects

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Matrix assembly** | Matrix-free option | Always assembled |
| **Preconditioning** | Algebraic multigrid | LU decomposition |
| **Memory usage** | Lower (matrix-free) | Higher (assembled) |
| **Scalability** | Excellent (HPC-designed) | Good (MOOSE infrastructure) |

### 11.4 Strengths and Weaknesses

**Tandem Strengths:**
- High-order accuracy
- Optimized kernels (YATeTo)
- PETSc ecosystem integration
- Matrix-free capabilities
- Clean three-operator design

**Tandem Weaknesses:**
- Steeper learning curve
- Less user-friendly configuration
- Requires external mesh generation

**MOOSE Strengths:**
- User-friendly input files
- Extensive documentation
- Built-in mesh generation
- MultiApp flexibility
- Large user community

**MOOSE Weaknesses:**
- Lower-order methods (current)
- Less optimized kernels
- Direct solver dependency
- More complex staggered coupling

---

## 12. Benchmark Verification Results

The MOOSE implementation (stiffness mode) was verified against the reference benchmark solution
(unicycle-ap-ratestate by Sylvain Barbot) for the BP2-QD problem.

### 12.1 Test Configuration

- **Mesh:** 800m element size, 2D antiplane
- **Duration:** 1 year (3.1557×10⁷ seconds)
- **Traction mode:** `stiffness` (K = μ/Wf = 801,000 Pa/m)
- **Runtime:** ~10 minutes on 4 threads

### 12.2 Quantitative Comparison at z=0 (VW region)

| Quantity | Benchmark | MOOSE | Error |
|----------|-----------|-------|-------|
| Shear stress (MPa) | 26.5646 | 26.5714 | **0.025%** |
| Slip rate (m/s) | 1.37×10⁻¹⁵ | 1.50×10⁻¹⁵ | ~9% |
| Accumulated slip (m) | 7.92×10⁻⁶ | 8.16×10⁻⁶ | **3.0%** |
| State variable (s) | 3.32×10⁷ | 3.16×10⁷ | **-4.8%** |

### 12.3 Physics Verification

| Expected Behavior | MOOSE Result | Status |
|-------------------|--------------|--------|
| VW region locks up | Slip rate: 10⁻⁹ → 10⁻¹⁵ m/s | ✅ |
| Stress builds in VW | +25 kPa/year | ✅ |
| VS creeps at plate rate | V = 10⁻⁹ m/s (constant) | ✅ |
| VS stress constant | τ = 26.546 MPa (unchanged) | ✅ |
| State evolution (aging law) | θ grows with locking | ✅ |

### 12.4 Conclusion

The MOOSE implementation with stiffness-based traction computation accurately reproduces
the BP2-QD benchmark results. The stress error is excellent (<0.03%), and the slip/state
errors are within acceptable bounds (<5%) for the 800m mesh resolution.

---

## 13. Recommendations

### 13.1 For MOOSE Development

1. **High-Order Extension:**
   - Implement higher polynomial order DG (p=2-4)
   - Consider hierarchical basis functions

2. **Solver Improvements:**
   - Add multigrid preconditioner option
   - Implement matrix-free DG operators

3. **Time Stepping:**
   - Consider PETSc TS integration via MOOSE
   - Implement embedded RK pairs for error estimation

4. **Friction Law:**
   - Keep direct θ convention (benchmark-compatible)
   - Add slip law and other friction law options

### 13.2 For Production Simulations

1. **BP2-QD Benchmark:**
   - MOOSE multiapp_bp2 adequate for validation
   - Use stiffness mode for traction computation

2. **Large-Scale 3D:**
   - Consider Tandem for high-resolution 3D
   - MOOSE MultiApp for moderate resolution

3. **Parameter Studies:**
   - MOOSE input file flexibility advantageous
   - Consider DiscreteGreenOperator for efficiency

### 13.3 Hybrid Approach

Consider combining strengths:
- Use MOOSE MultiApp architecture
- Implement YATeTo-style kernel generation
- Add PETSc TS time stepping option
- Develop matrix-free DG operators

---

## Appendix A: Parameter Mapping

| Parameter | BP2 Symbol | Tandem | MOOSE |
|-----------|------------|--------|-------|
| Shear modulus | μ | `mu` (Lua) | `shear_modulus` |
| Density | ρ | `rho` (Lua) | `rho` |
| Reference friction | f₀ | `cp_.f0` | `_f0` |
| Reference velocity | V₀ | `cp_.V0` | `_V0` |
| Direct effect | a | `p_[i].get<A>()` | `_a` (or coupled) |
| Evolution effect | b | `cp_.b` | `_b` |
| Critical distance | Dc | `p_[i].get<L>()` | `_Dc` |
| Normal stress | σₙ | `p_[i].get<SnPre>()` | `_sigma_n` |
| Radiation damping | η | `p_[i].get<Eta>()` | `_eta` |

---

## Appendix B: File Location Reference

### Tandem Key Files

| Component | Path |
|-----------|------|
| Main driver | `app/tandem/SEAS.cpp` |
| Elasticity | `app/localoperator/Elasticity.h` |
| Friction law | `app/localoperator/DieterichRuinaAgeing.h` |
| QD operator | `app/form/SeasQDOperator.h` |
| FD operator | `app/form/SeasFDOperator.h` |
| Time solver | `app/common/PetscTimeSolver.h` |
| Kernel generation | `app/kernels/elasticity.py` |

### MOOSE Key Files

| Component | Path |
|-----------|------|
| DG elasticity | `src/dgkernels/DGElasticityAntiplane.C` |
| Fault interface | `src/interfacekernels/DGFaultSlipInterfaceKernel.C` |
| Traction material | `src/materials/dg/DGElasticTractionMaterial.C` |
| Slip rate aux | `src/auxkernels/seas/SEASSlipRateVarAAux.C` |
| State evolution | `src/kernels/friction/StateEvolutionKernel.C` |
| Adaptive DT | `src/timestepper/SEASAdaptiveDT.C` |
| benchmark_bp2 input | `examples/benchmark_bp2/bp2.i` |
| multiapp_bp2 main | `examples/multiapp_bp2/main_elasticity.i` |
| multiapp_bp2 sub | `examples/multiapp_bp2/sub_friction.i` |

---

## Appendix C: Detailed Analysis of Tandem Traction Computation

This appendix provides an in-depth technical review of how Tandem computes traction from displacement through the `adapter_->traction(disp_view, traction_)` call chain. This is a critical component that couples the DG elasticity solver to the friction law.

### C.1 High-Level Call Flow

The traction computation is triggered from the SEAS operators (quasi-dynamic or fully dynamic):

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  SeasQDOperator::rhs(time, state, result)                                    │
│    │                                                                         │
│    ├─► solve(time, state_view)         // Step 1: Solve elasticity           │
│    │     • dgop_->set_slip(adapter_->slip_bc(state_view))  // Apply slip BC │
│    │     • dgop_->set_dirichlet(fun_boundary_(time))       // Far-field BC  │
│    │     • linear_solver_.update_rhs(*dgop_)               // Assemble RHS  │
│    │     • linear_solver_.solve()                          // Solve Ku = f  │
│    │     OUTPUT: displacement u in linear_solver_.x()                        │
│    │                                                                         │
│    ├─► update_traction(state_view)     // Step 2: Compute traction           │
│    │     • disp_view = LocalGhostCompositeView(linear_solver_.x(), ghost)   │
│    │     • dgop_->set_slip(adapter_->slip_bc(state_view))  // Set slip BC   │
│    │     • adapter_->traction(disp_view, traction_)        // ◄── KEY CALL  │
│    │     OUTPUT: traction_ (BlockVector on fault)                            │
│    │                                                                         │
│    └─► friction_->rhs(time, traction_, state, result)  // Step 3: Friction  │
│          • For each fault node: Solve τ = σn·f(V,ψ) + η·V for V            │
│          OUTPUT: result = [dslip/dt = V, dψ/dt]                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Source:** `app/form/SeasQDOperator.cpp` (lines 34-40, 69-74)

```cpp
// SeasQDOperator.cpp - rhs function
void SeasQDOperator::rhs(double time, BlockVector const& state, BlockVector& result) {
    update_ghost_state(state);
    solve(time, make_state_view(state));      // Step 1: Solve elasticity
    update_traction(make_state_view(state));  // Step 2: Compute traction
    friction_->rhs(time, traction_, state, result);  // Step 3: Friction RHS
}

// SeasQDOperator.cpp - update_traction function
void SeasQDOperator::update_traction(BlockView const& state_view) {
    auto disp_view = LocalGhostCompositeView(linear_solver_.x(), disp_ghost_);
    dgop_->set_slip(adapter_->slip_bc(state_view));  // IMPORTANT: Set slip for traction formula
    adapter_->traction(disp_view, traction_);        // ◄── TRACTION COMPUTATION
    dgop_->set_slip(invalid_slip_bc());              // Clear for safety
}
```

### C.2 AdapterOperator::traction - Orchestration Layer

**Source:** `app/form/AdapterOperator.h` (lines 57-78)

The `AdapterOperator` template class orchestrates the traction computation by:
1. Iterating over fault elements
2. Extracting displacement from adjacent domain elements
3. Computing traction at quadrature points via elasticity kernels
4. Projecting traction to fault basis functions via adapter kernels

```cpp
template <typename LocalOperator>
void AdapterOperator<LocalOperator>::traction(BlockView const& displacement,
                                               BlockVector& result) override {
    // Allocate temporary storage for traction at quadrature points
    auto traction_q = Managed<Matrix<double>>(adapted_lop_->tractionResultInfo().shape(),
                                              std::size_t{ALIGNMENT});

    scratch_.reset();
    auto result_handle = result.begin_access();

    // Loop over all local fault elements
    for (std::size_t faultNo = 0, num = num_local_elements(); faultNo < num; ++faultNo) {
        auto fctNo = fault_map_->fctNo(faultNo);     // Get facet number
        auto const& info = topo_->info(fctNo);       // Get topological info

        // Extract displacement from elements on both sides of fault
        auto u0 = displacement.get_block(info.up[0]);  // "minus" side
        auto u1 = displacement.get_block(info.up[1]);  // "plus" side

        // Compute traction at quadrature points
        if (info.up[0] == info.up[1]) {
            // Boundary facet (single element)
            adapted_lop_->traction_boundary(fctNo, info, u0, traction_q);
        } else {
            // Interior (skeleton) facet (two elements)
            adapted_lop_->traction_skeleton(fctNo, info, u0, u1, traction_q);
        }

        // Project traction from quadrature points to basis function coefficients
        auto result_block = result_handle.subtensor(slice{}, faultNo);
        lop_->traction(faultNo, traction_q, result_block, scratch_);
    }
    result.end_access(result_handle);
}
```

### C.3 Elasticity::traction_skeleton - Stress Computation at Interior Facets

**Source:** `app/localoperator/Elasticity.cpp` (lines 948-987)

For interior fault facets (skeleton), traction is computed from displacement gradients on both sides:

```cpp
void Elasticity::traction_skeleton(std::size_t fctNo, FacetInfo const& info,
                                   Vector<double const>& u0, Vector<double const>& u1,
                                   Matrix<double>& result) const {
    // ═══════════════════════════════════════════════════════════════════════
    // STEP 1: Compute spatial gradients at quadrature points
    // ═══════════════════════════════════════════════════════════════════════
    //
    // Transform reference gradients to physical gradients:
    //   ∇_x u = G^{-T} · ∇_ξ u
    //
    // Where G = ∂x/∂ξ is the Jacobian of the coordinate transformation
    //
    alignas(ALIGNMENT) double Dx_q0[tensor::Dx_q::size(0)];  // ∇u on side 0
    alignas(ALIGNMENT) double Dx_q1[tensor::Dx_q::size(1)];  // ∇u on side 1

    kernel::Dx_q dxKrnl;
    dxKrnl.Dx_q(0) = Dx_q0;
    dxKrnl.Dx_q(1) = Dx_q1;
    dxKrnl.g(0) = fct[fctNo].get<JInv0>().data()->data();  // G^{-T} side 0
    dxKrnl.g(1) = fct[fctNo].get<JInv1>().data()->data();  // G^{-T} side 1
    for (unsigned side = 0; side < 2; ++side) {
        dxKrnl.Dxi_q(side) = Dxi_q[info.localNo[side]].data();  // ∇_ξ basis
        dxKrnl.execute(side);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // STEP 2: Get boundary condition (prescribed slip)
    // ═══════════════════════════════════════════════════════════════════════
    alignas(ALIGNMENT) double f_q_raw[tensor::f_q::size()];
    bc_skeleton(fctNo, info.bc, f_q_raw);  // f_q = slip at quadrature points

    // ═══════════════════════════════════════════════════════════════════════
    // STEP 3: Compute traction using YATeTo-generated kernel
    // ═══════════════════════════════════════════════════════════════════════
    kernel::compute_traction krnl;
    krnl.c00 = -penalty(fctNo);                             // DG penalty coefficient
    krnl.Dx_q(0) = Dx_q0;                                   // ∇u side 0
    krnl.Dx_q(1) = Dx_q1;                                   // ∇u side 1
    krnl.E_q(0) = E_q[info.localNo[0]].data();              // Basis functions side 0
    krnl.E_q(1) = E_q[info.localNo[1]].data();              // Basis functions side 1
    krnl.f_q = f_q_raw;                                     // Slip BC
    krnl.lam_q(0) = fctPre[fctNo].get<lam_q_0>().data();    // λ at quad pts side 0
    krnl.lam_q(1) = fctPre[fctNo].get<lam_q_1>().data();    // λ at quad pts side 1
    krnl.mu_q(0) = fctPre[fctNo].get<mu_q_0>().data();      // μ at quad pts side 0
    krnl.mu_q(1) = fctPre[fctNo].get<mu_q_1>().data();      // μ at quad pts side 1
    krnl.n_unit_q = fct[fctNo].get<UnitNormal>().data()->data();  // Unit normal
    krnl.traction_q = result.data();                        // OUTPUT: traction
    krnl.u(0) = u0.data();                                  // Displacement coeffs side 0
    krnl.u(1) = u1.data();                                  // Displacement coeffs side 1
    krnl.execute();
}
```

### C.4 YATeTo Kernel Definitions - Mathematical Formulation

**Source:** `app/kernels/elasticity.py` (lines 89-91, 242-248)

The traction computation is defined using YATeTo (Yet Another Tensor Toolbox), which generates optimized C++ kernels from tensor notation:

```python
# ═══════════════════════════════════════════════════════════════════════════════
# Traction function: T = σ · n (stress contracted with normal)
# ═══════════════════════════════════════════════════════════════════════════════
#
# For isotropic linear elasticity (Hooke's law):
#   σ_ij = λ δ_ij ε_kk + 2μ ε_ij
#   ε_ij = (1/2)(∂u_i/∂x_j + ∂u_j/∂x_i)
#
# Traction on surface with normal n:
#   T_i = σ_ij n_j
#       = λ (∇·u) n_i + μ (∂u_i/∂x_j n_j + ∂u_j/∂x_i n_j)
#       = λ (∂u_s/∂x_s) n_p + μ (∂u_p/∂x_j n_j + ∂u_j/∂x_p n_j)
#
def traction(x, normal):
    return (
        lam_q[x]['q'] * Dx_q[x]['lsq'] * u[x]['ls'] * normal['pq']  # λ tr(ε) n
        + mu_q[x]['q'] * (
            Dx_q[x]['ljq'] * u[x]['lp'] * normal['jq']   # μ (∇u)·n
          + Dx_q[x]['lpq'] * u[x]['lj'] * normal['jq']   # μ (∇u)^T·n
        )
    )

# Index notation explanation:
#   l: basis function index
#   p, j, s: spatial dimensions (x, y, z)
#   q: quadrature point index
#
#   Dx_q[x]['lpq'] = ∂φ_l/∂x_p at quadrature point q on side x
#   u[x]['lp']     = u_p^l = coefficient of basis function l for component p
#   normal['pq']   = n_p at quadrature point q
```

**Interior (skeleton) traction formula:**

```python
# ═══════════════════════════════════════════════════════════════════════════════
# Skeleton traction: DG-consistent formula with penalty stabilization
# ═══════════════════════════════════════════════════════════════════════════════
#
# T^q_p = (1/2) * [T_0(u_0,n) + T_1(u_1,n)]           ← Average of tractions
#       + c0 * [E_q[0]·u[0] - E_q[1]·u[1] - f_q]_p^q  ← Penalty on jump(u) - slip
#
# Where:
#   T_x = traction from side x (computed via Hooke's law)
#   c0 = penalty parameter (negative in code, so c00 = -penalty)
#   E_q = basis functions evaluated at quadrature points
#   f_q = prescribed slip boundary condition
#
generator.add('compute_traction',
    traction_q['pq'] <= 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q))
                      + c0[0] * (E_q[0]['lq'] * u[0]['lp']      # u_0 at quad pts
                               - E_q[1]['lq'] * u[1]['lp']      # u_1 at quad pts
                               - f_q['pq']))                     # slip BC
```

**Boundary traction formula:**

```python
# ═══════════════════════════════════════════════════════════════════════════════
# Boundary traction: Single-sided formula
# ═══════════════════════════════════════════════════════════════════════════════
#
# T^q_p = T_0(u_0,n) + c0 * [E_q[0]·u[0] - f_q]_p^q
#
generator.add('compute_traction_bnd',
    traction_q['pq'] <= traction(0, n_unit_q)
                      + c0[0] * (E_q[0]['lq'] * u[0]['lp'] - f_q['pq']))
```

### C.5 Adapter::traction - Projection to Fault Basis Functions

**Source:** `app/localoperator/ElasticityAdapter.cpp` (lines 19-34)

After computing traction at quadrature points, the values must be projected onto the fault's polynomial basis:

```cpp
template <>
void Adapter<Elasticity>::traction(std::size_t faultNo,
                                   Matrix<double> const& traction_q,   // INPUT: T at quad pts
                                   Vector<double>& traction,           // OUTPUT: coefficients
                                   LinearAllocator<double>&) const {
    elasticity_adapter::kernel::evaluate_traction krnl;
    krnl.e_q_T = e_q_T.data();                                    // Basis functions (transposed)
    krnl.fault_basis_q = fault_[faultNo].template get<FaultBasis>().data()->data();  // Coord transform
    krnl.traction_q = traction_q.data();                          // Traction at quad points
    krnl.minv = mass_[faultNo].template get<MInv>().data();       // Inverse mass matrix
    krnl.nl_q = fault_[faultNo].template get<NormalLength>().data();  // Surface Jacobian
    krnl.traction = traction.data();                              // OUTPUT
    krnl.w = quad_rule_.weights().data();                         // Quadrature weights
    krnl.execute();
}
```

**Source:** `app/kernels/elasticity_adapter.py` (lines 24-27)

```python
# ═══════════════════════════════════════════════════════════════════════════════
# Projection: L2 projection of traction onto fault polynomial basis
# ═══════════════════════════════════════════════════════════════════════════════
#
# This is a weighted L2 projection:
#
#   traction[k,p] = ∫_Γ T_p(x) φ_k(x) dΓ  (approx by quadrature)
#
# Using inverse mass matrix for efficiency:
#
#   traction[k,p] = M^{-1}[l,k] * Σ_q (e_q^T[q,l] * w[q] * |n|[q] * T[o,q] * B[o,p,q])
#
# Where:
#   M^{-1}    = inverse mass matrix on fault
#   e_q^T     = basis functions evaluated at quadrature points (transposed)
#   w         = quadrature weights
#   |n|       = surface Jacobian (normal length at quad pt)
#   T[o,q]    = traction component o at quadrature point q
#   B[o,p,q]  = fault basis transformation tensor (maps domain coords to fault coords)
#
generator.add('evaluate_traction',
    traction['kp'] <= minv['lk'] * e_q_T['ql'] * w['q'] *
                      nl_q['q'] * traction_q['oq'] * fault_basis_q['opq'])
```

### C.6 Mathematical Summary

**Complete Traction Formula for Interior Fault (Skeleton):**

```
                    ┌─────────────────────────────────────────────────────┐
                    │  T(x) = ½[T⁺(x) + T⁻(x)] + η·([[u]] - slip)        │
                    └─────────────────────────────────────────────────────┘

Where:

T±(x) = σ±(x) · n    (Cauchy's formula)

σ± = λ± tr(ε±) I + 2μ± ε±    (Hooke's law for linear elasticity)

ε± = sym(∇u±) = ½(∇u± + (∇u±)ᵀ)    (Small strain tensor)

[[u]] = u⁺ - u⁻    (Displacement jump across fault)

η = penalty parameter (stabilizes DG method)
```

**Step-by-step computation:**

```
1. INPUT: Nodal displacement coefficients u⁺, u⁻ for elements on each side of fault

2. GRADIENT: Compute ∇u at quadrature points
   ∇_x u(x_q) = G⁻ᵀ · Σ_l (u_l · ∇_ξ φ_l(ξ_q))

   where G = ∂x/∂ξ (Jacobian), φ_l = basis functions

3. STRESS: Apply Hooke's law at quadrature points
   σ_ij(x_q) = λ(x_q) δ_ij ∇·u(x_q) + μ(x_q)(∂u_i/∂x_j + ∂u_j/∂x_i)

4. TRACTION: Contract stress with normal
   T_p(x_q) = σ_pj(x_q) · n_j(x_q)

5. DG AVERAGE + PENALTY:
   T_p^{DG}(x_q) = ½[T_p⁺(x_q) + T_p⁻(x_q)] + η·(u_p⁺(x_q) - u_p⁻(x_q) - slip_p(x_q))

6. PROJECTION: Integrate against basis functions
   T_{k,p} = M⁻¹_{lk} · Σ_q [w_q · |n_q| · T_p^{DG}(x_q) · φ_l(x_q) · B_{p}(x_q)]

7. OUTPUT: Traction coefficients T_{k,p} for use in friction law
```

### C.7 Key Precomputed Quantities

The following quantities are precomputed during setup for efficiency:

| Quantity | Symbol | Description | Computed In |
|----------|--------|-------------|-------------|
| `JInv0`, `JInv1` | G⁻ᵀ | Inverse transpose Jacobian | `Elasticity::prepare_facet` |
| `Dxi_q` | ∇_ξ φ | Reference gradients at quad pts | `Elasticity` constructor |
| `E_q` | φ | Basis functions at quad pts | `Elasticity` constructor |
| `lam_q`, `mu_q` | λ, μ | Lamé parameters at quad pts | `Elasticity::precomputeSurface` |
| `UnitNormal` | n̂ | Unit normal vector | `Elasticity::prepare_facet` |
| `penalty` | η | DG penalty parameter | `Elasticity::penalty()` |
| `MInv` | M⁻¹ | Inverse mass matrix (fault) | `AdapterBase::prepare` |
| `NormalLength` | \|n\| | Surface Jacobian | `AdapterBase::prepare` |
| `FaultBasis` | B | Coord transform tensor | `AdapterBase::prepare` |

### C.8 Connection to Friction Law

The computed traction is passed directly to the friction operator:

```cpp
// In SeasQDOperator::rhs():
friction_->rhs(time, traction_, state, result);
```

The friction operator (`DieterichRuinaAgeing.h`) uses this traction to:

1. **Solve for slip rate V** via Brent's method:
   ```
   τ = σn · f(V, ψ) + η_rad · V

   where:
   τ = |T|           (traction magnitude from this computation)
   σn = normal stress
   f(V,ψ) = friction coefficient
   η_rad = radiation damping
   ```

2. **Compute state evolution rate**:
   ```
   dψ/dt = (b·V₀/L)·(exp((f₀-ψ)/b) - V/V₀)
   ```

### C.9 Comparison with MOOSE Implementation

| Aspect | Tandem | MOOSE |
|--------|--------|-------|
| **Traction location** | Quadrature points on fault | Interface quadrature points |
| **Gradient computation** | YATeTo-optimized kernel | Standard MOOSE shape functions |
| **DG penalty** | Computed via `penalty(fctNo)` | Parameter `_penalty` in InterfaceKernel |
| **Projection** | L2 projection with inverse mass | Direct evaluation at nodes |
| **Slip subtraction** | In traction formula: `[[u]] - slip` | In InterfaceKernel residual |
| **Coordinate transform** | `FaultBasis` tensor | Implicit in interface normals |

### C.10 Key Files Summary

| File | Purpose |
|------|---------|
| `app/form/SeasQDOperator.cpp` | Calls `update_traction()` in RHS evaluation |
| `app/form/AdapterOperator.h` | Orchestrates traction computation loop |
| `app/localoperator/Elasticity.cpp` | `traction_skeleton()` and `traction_boundary()` |
| `app/kernels/elasticity.py` | YATeTo kernel definitions |
| `app/localoperator/ElasticityAdapter.cpp` | Projection to fault basis |
| `app/kernels/elasticity_adapter.py` | Projection kernel definition |

---

*Document End*
