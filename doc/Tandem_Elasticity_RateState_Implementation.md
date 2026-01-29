# Tandem: Elasticity and Rate-and-State Friction Implementation

This document describes how the Tandem code implements elasticity coupled with rate-and-state friction for earthquake simulations (SEAS - Sequences of Earthquakes and Aseismic Slip).

## 1. Overview

**Tandem** is a high-order Discontinuous Galerkin (DG) finite element code for earthquake and aseismic slip simulation. It uses:
- Symmetric Interior Penalty Galerkin (SIPG) method for elasticity
- PETSc for time stepping and linear solvers
- YATeTo for optimized tensor kernel generation

**Key Dependencies:**
- PETSc (linear/nonlinear solvers, time stepping)
- Eigen (linear algebra)
- YATeTo (tensor kernel code generation)
- METIS/ParMETIS (mesh partitioning)
- MPI (distributed memory parallelization)

---

## 2. Elasticity Implementation

### 2.1 Discretization Method

Tandem uses the **Symmetric Interior Penalty Galerkin (SIPG)** method for discretizing the elastodynamic equations on unstructured simplicial meshes.

**Location:** `/app/localoperator/Elasticity.h`

**Governing Equations:**
```
ρ ∂²u/∂t² = ∇·σ + f_body

σ = λ tr(ε(u)) I + 2μ ε(u)

ε(u) = ½(∇u + (∇u)ᵀ)
```

### 2.2 Key Class Structure

```cpp
class Elasticity : public DGCurvilinearCommon<DomainDimension>
```

**Precomputed Matrices (reference element):**
- `MhatInv`: Inverse mass matrix
- `E_Q`, `E_Q_T`: Basis function evaluations at quadrature points
- `Dxi_Q`, `Dxi_Q_120`: Basis gradients at quadrature points

**Material Parameters:**
- λ (Lamé's first parameter)
- μ (shear modulus)
- ρ (density)

### 2.3 Assembly Methods

| Method | Description |
|--------|-------------|
| `prepare_volume()` | Precompute volume contributions |
| `prepare_skeleton()` | Precompute interior facet contributions |
| `prepare_boundary()` | Precompute boundary facet contributions |
| `assemble_volume()` | Assemble volume stiffness |
| `assemble_skeleton()` | Assemble interior penalty terms |
| `assemble_boundary()` | Assemble boundary conditions |
| `rhs_volume()` | Volume RHS contribution |
| `rhs_skeleton()` | Interior facet RHS |
| `rhs_boundary()` | Boundary RHS |

### 2.4 Kernel Generation

Tandem uses **YATeTo** (Yet Another Tensor Toolbox) to auto-generate optimized tensor kernels.

**Location:** `/app/kernels/elasticity.py`

Generated kernels include:
- `volumeOp`: Volume contribution with λ and μ
- `surfaceOp`, `surfaceOpBnd`: Facet contributions
- `lift_ip`, `lift_skeleton`, `lift_boundary`: Interior penalty lifting

---

## 3. Rate-and-State Friction Implementation

### 3.1 Class Hierarchy

```
RateAndStateBase (base class)
    └── RateAndState<Law> (templated wrapper)
            └── DieterichRuinaAgeing (specific friction law)
```

**Locations:**
- `/app/localoperator/RateAndStateBase.h`
- `/app/localoperator/RateAndState.h`
- `/app/localoperator/DieterichRuinaAgeing.h`

### 3.2 Dieterich-Ruina Aging Law

**State Evolution Equation:**
```
dψ/dt = (b·V₀/L) · (exp((f₀-ψ)/b) - V/V₀)
```

**Friction Coefficient:**
```
f(V, ψ) = a · sinh⁻¹((V/(2V₀)) · exp(ψ/a))
```

**Shear Stress Balance:**
```
τ = σₙ · f(V, ψ) + η·V
```

Where:
- ψ: State variable
- V: Slip rate magnitude
- V₀: Reference slip rate
- L: State evolution distance (Dc)
- a, b: Rate-and-state friction parameters
- f₀: Reference friction coefficient
- η: Radiation damping coefficient
- σₙ: Normal stress (positive in compression)
- τ: Shear stress magnitude

### 3.3 Friction Parameters

```cpp
struct ConstantParams {
    double V0;     // Reference slip rate
    double b;      // Friction evolution parameter
    double f0;     // Reference friction coefficient
};

struct Params {
    double a;      // Direct effect parameter
    double eta;    // Radiation damping coefficient
    double L;      // State evolution distance
    double sn_pre; // Pre-stress (normal)
    std::array<double, TangentialComponents> tau_pre;  // Pre-stress (shear)
    std::array<double, TangentialComponents> Vinit;    // Initial slip rate
    std::array<double, TangentialComponents> Sinit;    // Initial slip
};
```

### 3.4 Slip Rate Computation

The slip rate V is found by solving the nonlinear equation:
```
|τ̂| = σ̂ₙ · f(|V|, ψ) + η·|V|
```

Where:
- τ̂ = τ + τ_pre + η·V (absolute shear traction)
- σ̂ₙ = -σₙ + σₙ_pre (absolute normal stress)

**Solution Method:** Brent's method root finding (`zeroIn()` function)

---

## 4. Elasticity-Friction Coupling

### 4.1 Three-Operator Architecture

Tandem uses a three-operator system for coupling:

```
┌─────────────────────────────────────────────────────────┐
│                    SEAS Operator                         │
│  ┌─────────────┐  ┌──────────────┐  ┌────────────────┐  │
│  │   Domain    │  │   Adapter    │  │   Friction     │  │
│  │  Operator   │◄─┤   Operator   │◄─┤   Operator     │  │
│  │ (Elasticity)│  │  (Coupling)  │  │ (Rate-State)   │  │
│  └─────────────┘  └──────────────┘  └────────────────┘  │
└─────────────────────────────────────────────────────────┘
```

1. **Domain Operator** (`SeasFDOperator`): Solves elastodynamics
2. **Adapter Operator**: Computes traction from displacement, applies slip BC
3. **Friction Operator**: Evolves friction state and computes slip rate

### 4.2 Coupling Workflow

```
State Vector = [v, u, s]  (velocity, displacement, fault state)

For each time step:
    1. Compute displacement from velocity: du/dt = v

    2. Compute velocity from elasticity:
       dv/dt = (1/ρ)(∇·σ + f)  with slip boundary condition

    3. Compute traction on fault:
       t = σ·n  from displacement gradient

    4. Evolve fault state:
       ds/dt = friction.rhs(traction, state)
       └─ Extract τ, σₙ from traction
       └─ Solve for V(τ, σₙ, ψ)
       └─ Compute dψ/dt
```

### 4.3 Adapter Implementation

**Location:** `/app/localoperator/Adapter.h`, `ElasticityAdapter.cpp`

Key methods:
- `traction()`: Computes traction from stress tensor at fault
- `slip()`: Maps slip displacement to boundary condition

The adapter:
1. Computes stress σ from displacement gradient
2. Projects traction t = σ·n onto fault surface
3. Transforms to fault-local coordinates (strike, dip)
4. Uses L2 projection with fault mass matrix

---

## 5. Solvers

### 5.1 Time Integration

**Location:** `/app/common/PetscTimeSolver.h`

Tandem wraps PETSc's TS (Time Stepping) module:

```cpp
template <std::size_t NumStateVecs>
class PetscTimeSolver
```

**Available Integrators (via PETSc):**
- Explicit Runge-Kutta: RK2, RK3, RK4
- Implicit: BDF, Theta methods
- Adaptive time stepping with error control

**Key Methods:**
- `TSSetRHSFunction()`: Sets the right-hand side function
- `solve()`: Advances solution to specified time
- `set_monitor()`: Sets output monitoring callback

**CFL Stability:**
For explicit schemes, the time step is limited by:
```
Δt ≤ CFL × h / c_max
```
Where c_max is the maximum wave speed (P-wave or S-wave).

### 5.2 Linear Solver

**Location:** `/app/common/PetscLinearSolver.h`

```cpp
class PetscLinearSolver
```

**Features:**
- Matrix-free and matrix-based modes
- Multigrid preconditioners
- Iterative Krylov methods (CG, GMRES, MINRES)

**Multigrid Configuration (`MGConfig.h`):**
- Coarse level polynomial degree
- Strategies: TwoLevel, Logarithmic, Full

### 5.3 Nonlinear Solver for Friction

The slip rate V is found by root finding:

```cpp
// From DieterichRuinaAgeing.h
double g = zeroIn(VLo, VHi, [&](double V) {
    return tau_abs - sn_hat * f(V, psi) - eta * V;
}, 1.0e-14);
```

Uses **Brent's method** for robust root finding within bounds [VLo, VHi].

---

## 6. Simulation Modes

Tandem supports three SEAS simulation modes:

### 6.1 Fully Dynamic (FD)
- Time-dependent elastodynamics
- Explicit time stepping with CFL constraint
- Full inertial effects included
- State: [velocity, displacement, fault state]

### 6.2 Quasi-Dynamic (QD)
- Static elasticity at each time step
- Radiation damping approximates inertia
- Implicit integration of friction
- Faster than FD for slow slip

### 6.3 Quasi-Dynamic with Discrete Green's Functions (QDGreen)
- Pre-computed Green's function for repeated scenarios
- Fastest for parameter studies
- Checkpointing support

---

## 7. Configuration

### 7.1 TOML Parameters

```toml
final_time = 0.01
resolution = 1
mode = "FD"              # FullyDynamic, QuasiDynamic, QuasiDynamicDiscreteGreen
type = "elasticity"
ref_normal = [1, 0, 0]   # Fault normal direction
cfl = 1.0                # CFL parameter
matrix_free = false      # Matrix-free operators
mg_strategy = "TwoLevel" # Multigrid strategy
```

### 7.2 Lua Material/Friction Definition

```lua
function scenario:rho(x, y, z)    -- Density
function scenario:mu(x, y, z)     -- Shear modulus
function scenario:lam(x, y, z)    -- Lamé parameter
function scenario:a(x, y, z)      -- Rate-state 'a'
function scenario:b(x, y, z)      -- Rate-state 'b'
function scenario:L(x, y, z)      -- State evolution distance
function scenario:sn_pre(x, y, z) -- Initial normal stress
```

---

## 8. Key File Locations

| Component | File Path |
|-----------|-----------|
| Elasticity solver | `app/localoperator/Elasticity.h` |
| Friction law | `app/localoperator/DieterichRuinaAgeing.h` |
| Rate-state wrapper | `app/localoperator/RateAndState.h` |
| Elastic-friction coupling | `app/localoperator/Adapter.h` |
| Fully-dynamic operator | `app/form/SeasFDOperator.h/.cpp` |
| Friction operator | `app/form/FrictionOperator.h` |
| Time solver | `app/common/PetscTimeSolver.h` |
| Linear solver | `app/common/PetscLinearSolver.h` |
| Main SEAS driver | `app/tandem/SEAS.cpp` |
| Kernel generation | `app/kernels/elasticity.py` |

---

## 9. Summary

**Elasticity:**
- SIPG DG method on unstructured simplicial meshes
- Matrix-free option with YATeTo optimized kernels
- Multigrid preconditioning for linear solves

**Rate-and-State Friction:**
- Dieterich-Ruina aging law
- State variable ψ evolved with slip rate
- Brent's method for slip rate inversion

**Solvers:**
- PETSc TS for time integration (explicit RK or implicit)
- PETSc KSP for linear systems (CG, GMRES with multigrid)
- Brent's method for nonlinear friction solve

**Coupling:**
- Three-operator architecture (Domain, Adapter, Friction)
- Traction computed from displacement at fault
- Slip rate applied as boundary condition to elasticity
