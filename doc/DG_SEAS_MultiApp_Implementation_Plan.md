# MultiApp SEAS Implementation Plan for MOOSE

## Overview

This document outlines an updated implementation plan for SEAS (Sequences of Earthquakes and Aseismic Slip) simulations in MOOSE using a **MultiApp architecture** inspired by Tandem's three-operator design. This approach separates elasticity and friction into distinct applications with data transfer between them.

**Key Design Principle:** Like Tandem, we decouple the problem into:
1. **Elasticity App (MainApp)**: Solves elastodynamics/elastostatics with prescribed slip BC
2. **Friction App (SubApp)**: Evolves rate-and-state friction law, computes slip rate

---

## 1. Comparison: Tandem vs Current MOOSE vs Proposed MultiApp

### 1.1 Tandem Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                    SeasFDOperator (Main)                      │
│  ┌────────────────┐  ┌─────────────┐  ┌───────────────────┐  │
│  │ DGOperator     │  │  Adapter    │  │ FrictionOperator  │  │
│  │ (Elasticity)   │◄─┤  Operator   │◄─┤ (RateState)       │  │
│  │                │  │  (Coupling) │  │                   │  │
│  └────────────────┘  └─────────────┘  └───────────────────┘  │
└──────────────────────────────────────────────────────────────┘

Data Flow:
  1. Elasticity receives slip → computes displacement
  2. Adapter computes traction from displacement at fault
  3. Friction receives traction → solves for V, evolves θ → produces slip
```

### 1.2 Current Single-App MOOSE (Staggered)

```
┌─────────────────────────────────────────────────────────────┐
│                    Single MOOSE App                          │
│  ┌─────────────────┐   ┌─────────────────────────────────┐  │
│  │ DGKernels       │   │ AuxKernels (on fault boundary)  │  │
│  │ (Elasticity)    │   │ - SEASSlipAux                   │  │
│  │                 │   │ - SEASStateAux                  │  │
│  │ InterfaceKernel │◄──│ - SEASTractionAux              │  │
│  │ (FaultSlip)     │   │ - SEASSlipRateAux              │  │
│  └─────────────────┘   └─────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────┘

Limitation: Friction is solved on boundary of bulk mesh, not independent
```

### 1.3 Proposed MultiApp Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                         MOOSE MultiApp                               │
│                                                                      │
│  ┌──────────────────────────────────┐                               │
│  │        MainApp (Elasticity)       │                               │
│  │  - 2D/3D bulk mesh                │      ┌───────────────────┐   │
│  │  - DGKernels (elasticity)         │      │  Transfers        │   │
│  │  - Slip BC from SubApp            │◄─────┤  - Slip → BC      │   │
│  │  - Computes traction at fault     │─────►│  - Traction → Sub │   │
│  └──────────────────────────────────┘      └───────────────────┘   │
│                                                      │               │
│                                                      ▼               │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                    SubApp (Friction)                          │   │
│  │  - 1D/2D fault mesh (lower dimension)                         │   │
│  │  - Rate-State ODEs: dψ/dt, ds/dt                             │   │
│  │  - Receives traction from MainApp                             │   │
│  │  - Computes slip rate V, state θ                             │   │
│  │  - Sends slip back to MainApp                                 │   │
│  └──────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 2. Mathematical Formulation

### 2.1 Elasticity Problem (MainApp)

**Elastodynamics (Fully Dynamic):**
```
ρ ∂²u/∂t² = ∇·σ + f          in Ω
σ = C:ε(u)                    (constitutive)
u = g_D                       on Γ_D (Dirichlet)
σ·n = t                       on Γ_N (Neumann)
[[u]] = s_prescribed          on Γ_F (Fault - from SubApp)
```

**Quasi-Static (Quasi-Dynamic):**
```
∇·σ = 0                       in Ω
[[u]] = s_prescribed          on Γ_F (Fault - from SubApp)
```

The elasticity problem is **LINEAR** when slip is prescribed (from SubApp).

### 2.2 Friction Problem (SubApp)

**State Variables on Fault:**
- `s(t)`: Slip (displacement jump across fault)
- `V(t)`: Slip rate = ds/dt
- `ψ(t)` or `θ(t)`: State variable

**Governing Equations (Dieterich-Ruina Aging Law):**
```
Slip rate:     ds/dt = V

State evolution (aging law):
               dψ/dt = (b·V₀/L)·(exp((f₀-ψ)/b) - V/V₀)
         or    dθ/dt = 1 - V·θ/L

Traction balance:
               |τ̂| = σ̂ₙ·f(|V|, ψ) + η·|V|

where:
               τ̂ = τ_elastic + τ_pre         (from MainApp + initial stress)
               σ̂ₙ = |σₙ_pre|                 (effective normal stress)
               η = μ/(2·cs)                  (radiation damping)
```

**Regularized Friction Coefficient (Tandem Eq. 7):**
```
f(V, ψ) = a·sinh⁻¹((V/(2V₀))·exp(ψ/a))
```

### 2.3 Coupling: Adapter/Transfer

**MainApp → SubApp (Traction Transfer):**
```
τ_elastic = μ·{{∂u/∂n}}     (average of displacement gradient at fault)
```
Interpolate from bulk mesh boundary to fault mesh.

**SubApp → MainApp (Slip Transfer):**
```
s_prescribed = ∫V dt        (integrated slip from friction app)
```
Interpolate from fault mesh to bulk mesh boundary condition.

---

## 3. Implementation Architecture

### 3.1 File Structure

```
farms_rsf_explicit/
├── include/
│   ├── multiapps/
│   │   └── SEASFrictionApp.h           # SubApp declaration
│   ├── transfers/
│   │   ├── SEASTractionTransfer.h      # Main→Sub: traction transfer
│   │   └── SEASSlipTransfer.h          # Sub→Main: slip transfer
│   ├── kernels/
│   │   └── friction/
│   │       ├── SlipRateKernel.h        # ODE kernel for ds/dt = V
│   │       ├── StateEvolutionKernel.h  # ODE kernel for dθ/dt
│   │       └── TractionBalanceKernel.h # Nonlinear solve for V
│   ├── auxkernels/
│   │   └── friction/
│   │       ├── FrictionCoefficientAux.h  # Compute f(V,θ)
│   │       └── TractionAux.h             # Store received traction
│   ├── bcs/
│   │   └── PrescribedSlipBC.h          # Apply slip from SubApp
│   ├── userobjects/
│   │   └── FaultGeometry.h             # Fault coordinate system
│   └── timeintegrators/
│       └── SEASTimeIntegrator.h        # Custom integrator for friction
├── src/
│   └── [corresponding .C files]
└── examples/
    └── multiapp_bp2/
        ├── main_elasticity.i           # MainApp input
        ├── sub_friction.i              # SubApp input
        └── mesh/
            ├── bulk_mesh.msh           # 2D bulk mesh
            └── fault_mesh.msh          # 1D fault mesh
```

### 3.2 Class Hierarchy

```
MultiApp (MOOSE)
    └── TransientMultiApp
            └── SEASMultiApp (optional custom)

MultiAppTransfer (MOOSE)
    ├── MultiAppNearestNodeTransfer (for traction)
    ├── MultiAppProjectionTransfer (for slip)
    └── SEASTractionTransfer (custom if needed)

ScalarKernel (MOOSE)
    └── StateEvolutionKernel (NEW) - dθ/dt = g(V, θ)

AuxKernel (MOOSE)
    ├── SlipRateAux (solve V from τ balance)
    └── FrictionCoefficientAux (compute f(V,θ))

IntegratedBC (MOOSE)
    └── PrescribedSlipBC (apply [[u]] = s from SubApp)
```

---

## 4. Detailed Component Specifications

### 4.1 MainApp: Elasticity Application

**Input File Structure (`main_elasticity.i`):**

```
[Mesh]
  # 2D bulk mesh with fault boundary defined
  type = FileMesh
  file = bulk_mesh.msh
[]

[Variables]
  [u]  # Displacement (or w for antiplane)
    family = MONOMIAL  # For DG
    order = FIRST
  []
[]

[DGKernels]
  [dg_elasticity]
    type = DGElasticityAntiplane  # or DGElasticity3D
    variable = u
    epsilon = 1.0
    sigma = 6.0
  []
[]

[InterfaceKernels]
  [fault_slip]
    type = DGPrescribedSlipInterfaceKernel
    variable = u
    neighbor_var = u
    boundary = fault
    slip_variable = slip_from_subapp  # Transferred from SubApp
    penalty = 1e12
  []
[]

[AuxVariables]
  # Slip received from friction SubApp
  [slip_from_subapp]
    order = CONSTANT
    family = MONOMIAL
    boundary = fault
  []
  # Traction to send to SubApp
  [traction_to_subapp]
    order = CONSTANT
    family = MONOMIAL
    boundary = fault
  []
[]

[AuxKernels]
  [compute_traction]
    type = ElasticTractionAux
    variable = traction_to_subapp
    displacement = u
    shear_modulus = 32.04e9
    boundary = fault
    execute_on = 'TIMESTEP_END'
  []
[]

[MultiApps]
  [friction]
    type = TransientMultiApp
    input_files = sub_friction.i
    positions = '0 0 0'  # Origin of fault mesh
    execute_on = 'TIMESTEP_BEGIN'  # Friction updated before elasticity
  []
[]

[Transfers]
  # Send traction to SubApp
  [traction_to_friction]
    type = MultiAppNearestNodeTransfer
    to_multi_app = friction
    source_variable = traction_to_subapp
    variable = traction_from_main
    source_boundary = fault
  []
  # Receive slip from SubApp
  [slip_from_friction]
    type = MultiAppNearestNodeTransfer
    from_multi_app = friction
    source_variable = slip
    variable = slip_from_subapp
    source_boundary = fault_surface
  []
[]
```

### 4.2 SubApp: Friction Application

**Input File Structure (`sub_friction.i`):**

```
[Mesh]
  # 1D mesh representing fault (z-direction for BP2)
  type = GeneratedMesh
  dim = 1
  nx = 50              # 50 nodes along fault
  xmin = 0             # z = 0 (free surface)
  xmax = 40000         # z = 40 km (fault depth Wf)
[]

[Variables]
  # State variable θ (evolved via ODE)
  [state_variable]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  # Traction received from MainApp
  [traction_from_main]
    order = FIRST
    family = LAGRANGE
  []
  # Slip rate (solved algebraically)
  [slip_rate]
    order = FIRST
    family = LAGRANGE
  []
  # Accumulated slip
  [slip]
    order = FIRST
    family = LAGRANGE
  []
  # Friction coefficient
  [friction_coeff]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Time derivative of state variable
  [state_time]
    type = TimeDerivative
    variable = state_variable
  []
  # State evolution: dθ/dt = 1 - V*θ/Dc
  [state_evolution]
    type = StateEvolutionKernel
    variable = state_variable
    slip_rate = slip_rate
    Dc = 0.004
  []
[]

[AuxKernels]
  # Solve for slip rate from traction balance
  [solve_slip_rate]
    type = SlipRateAux
    variable = slip_rate
    traction = traction_from_main
    state_variable = state_variable
    a = 0.025
    b = 0.015
    Dc = 0.004
    f0 = 0.6
    V0 = 1e-6
    sigma_n = 50e6
    eta = 4.634e6  # μ/(2*cs)
    tau_pre = 26.5e6
    execute_on = 'TIMESTEP_BEGIN'
  []

  # Integrate slip: s += V * dt
  [integrate_slip]
    type = SlipIntegrationAux
    variable = slip
    slip_rate = slip_rate
    execute_on = 'TIMESTEP_END'
  []

  # Compute friction coefficient for output
  [compute_friction]
    type = FrictionCoefficientAux
    variable = friction_coeff
    slip_rate = slip_rate
    state_variable = state_variable
    a = 0.025
    b = 0.015
    f0 = 0.6
    V0 = 1e-6
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  [state_ic]
    type = ConstantIC
    variable = state_variable
    value = 4e6  # θ₀ = Dc/Vinit
  []
  [slip_ic]
    type = ConstantIC
    variable = slip
    value = 0
  []
  [slip_rate_ic]
    type = ConstantIC
    variable = slip_rate
    value = 1e-9  # Vinit
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  dt = 1e6
[]
```

### 4.3 Transfer Mechanisms

#### Option A: Use Built-in MOOSE Transfers

```
[Transfers]
  # Traction: MainApp boundary → SubApp volume
  [traction_transfer]
    type = MultiAppNearestNodeTransfer
    to_multi_app = friction
    source_variable = traction_to_subapp
    variable = traction_from_main
    source_boundary = fault
  []

  # Slip: SubApp volume → MainApp boundary
  [slip_transfer]
    type = MultiAppNearestNodeTransfer
    from_multi_app = friction
    source_variable = slip
    variable = slip_from_subapp
  []
[]
```

#### Option B: Custom SEASTransfer (for more control)

**SEASTractionTransfer** - Projects traction from 2D boundary to 1D mesh:
- Handles coordinate transformation (2D boundary coords → 1D fault coords)
- Averages traction from both sides of fault
- Handles element mismatch between meshes

**SEASSlipTransfer** - Projects slip from 1D mesh back to 2D boundary:
- Maps 1D fault coordinates to 2D boundary nodes
- Handles sign convention (slip direction)

### 4.4 New Kernels for SubApp

#### StateEvolutionKernel

**Purpose:** Implement dθ/dt = 1 - V·θ/Dc (aging law)

```cpp
class StateEvolutionKernel : public Kernel
{
public:
  static InputParameters validParams();
  StateEvolutionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableValue & _slip_rate;  // V (from AuxVariable)
  const Real _Dc;                     // Critical slip distance
};
```

**Residual:**
```
R = -1 + V*θ/Dc  (when combined with TimeDerivative gives dθ/dt = 1 - V*θ/Dc)
```

#### SlipRateAux

**Purpose:** Solve traction balance for V using Brent's method

```cpp
class SlipRateAux : public AuxKernel
{
protected:
  virtual Real computeValue() override;
  Real solveSlipRate(Real tau, Real theta) const;
  Real frictionCoefficient(Real V, Real theta) const;

private:
  const VariableValue & _traction;      // τ from MainApp
  const VariableValue & _state_variable; // θ
  // Rate-state parameters
  const Real _a, _b, _Dc, _f0, _V0;
  const Real _sigma_n;  // Normal stress
  const Real _eta;      // Radiation damping
  const Real _tau_pre;  // Initial shear stress
};
```

### 4.5 Prescribed Slip Interface Kernel

**DGPrescribedSlipInterfaceKernel:**
```cpp
class DGPrescribedSlipInterfaceKernel : public InterfaceKernel
{
protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

private:
  const VariableValue & _slip_prescribed;  // From SubApp
  const Real _penalty;
};
```

This enforces: `[[u]] = u_elem - u_neighbor = slip_prescribed`

---

## 5. Execution Flow

### 5.1 Time Step Execution Order

```
For each time step n → n+1:

┌─────────────────────────────────────────────────────────────┐
│ 1. TIMESTEP_BEGIN: SubApp Friction Update                    │
│    ├─ Receive τⁿ from MainApp (from previous step)          │
│    ├─ Solve for Vⁿ⁺¹ from: τⁿ = σₙ·f(Vⁿ⁺¹, θⁿ) + η·Vⁿ⁺¹    │
│    ├─ Update state: θⁿ⁺¹ = (θⁿ + Δt)/(1 + Vⁿ⁺¹·Δt/Dc)      │
│    └─ Update slip: sⁿ⁺¹ = sⁿ + Vⁿ⁺¹·Δt                      │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. TRANSFER: SubApp → MainApp                                │
│    └─ Send sⁿ⁺¹ to MainApp as slip BC                        │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. MAIN SOLVE: MainApp Elasticity                            │
│    ├─ Apply slip BC: [[u]] = sⁿ⁺¹                            │
│    ├─ Solve linear elasticity: ∇·σ = 0 (or ρü = ∇·σ)        │
│    └─ Compute uⁿ⁺¹                                           │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. TIMESTEP_END: Compute Traction                            │
│    └─ τⁿ⁺¹ = μ·{{∂uⁿ⁺¹/∂n}} at fault boundary               │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. TRANSFER: MainApp → SubApp                                │
│    └─ Send τⁿ⁺¹ to SubApp for next time step                 │
└─────────────────────────────────────────────────────────────┘
```

### 5.2 MultiApp Execution Configuration

```
[MultiApps]
  [friction]
    type = TransientMultiApp
    input_files = sub_friction.i
    execute_on = 'TIMESTEP_BEGIN'
    sub_cycling = false  # Keep same dt as main
  []
[]

[Transfers]
  [traction_to_sub]
    type = MultiAppNearestNodeTransfer
    to_multi_app = friction
    source_variable = traction
    variable = traction_from_main
    execute_on = 'TIMESTEP_BEGIN'  # Send before SubApp solves
  []
  [slip_from_sub]
    type = MultiAppNearestNodeTransfer
    from_multi_app = friction
    source_variable = slip
    variable = slip_from_subapp
    execute_on = 'TIMESTEP_BEGIN'  # Receive after SubApp solves
  []
[]
```

---

## 6. Comparison with Tandem Solvers

### 6.1 PETSc Time Integrator (Tandem)

Tandem uses PETSc TS with RK integrators:
```cpp
// Tandem: PetscTimeSolver.h
TSSetRHSFunction(ts, NULL, rhsFunction, ctx);
TSSolve(ts, solution);
```

**MOOSE Equivalent:**
- Use `Transient` Executioner with explicit time integrators
- Or implement custom `SEASTimeIntegrator` for specialized schemes

### 6.2 Linear Solver (Tandem)

Tandem uses PETSc KSP with multigrid:
```cpp
// Tandem: PetscLinearSolver.h
KSPSetOperators(ksp, A, P);
KSPSolve(ksp, b, x);
```

**MOOSE Equivalent:**
```
[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
```

### 6.3 Nonlinear Friction Solve (Tandem)

Tandem uses Brent's method for slip rate:
```cpp
// Tandem: DieterichRuinaAgeing.h
double V = zeroIn(VLo, VHi, [&](double V) {
    return tau - sn * f(V, psi) - eta * V;
}, tol);
```

**MOOSE Equivalent (already implemented):**
```cpp
// farms_rsf_explicit: BrentRootFinder.h
Real V = BrentRootFinder::zeroIn(VLo, VHi,
    [&](Real V) { return tau - sigma_n * f(V, theta) - eta * V; },
    tol);
```

---

## 7. Advantages of MultiApp Approach

| Aspect | Single-App (Current) | MultiApp (Proposed) |
|--------|---------------------|---------------------|
| **Mesh Flexibility** | Same mesh for bulk & fault | Different meshes, resolutions |
| **Code Modularity** | Friction embedded in AuxKernels | Separate friction physics |
| **Time Stepping** | Single dt for everything | Can sub-cycle friction |
| **Debugging** | Everything coupled | Can test apps independently |
| **Scalability** | All DOFs in one system | Distributed across apps |
| **Law Swapping** | Modify AuxKernels | Swap SubApp input file |
| **3D Extension** | 3D bulk + 2D fault boundary | 3D bulk + 2D fault mesh |

---

## 8. Implementation Phases

### Phase 1: Basic MultiApp Infrastructure
1. Create skeleton SubApp for friction (`sub_friction.i`)
2. Implement basic transfers (traction, slip)
3. Test with constant slip rate (bypass friction solve)
4. Verify traction computation matches single-app

### Phase 2: Friction SubApp Implementation
1. Implement `StateEvolutionKernel` for θ evolution
2. Implement `SlipRateAux` with Brent's method
3. Implement `SlipIntegrationAux` for s update
4. Test friction SubApp standalone with prescribed traction

### Phase 3: Coupling Verification
1. Connect MainApp and SubApp with transfers
2. Verify two-way data exchange
3. Test with BP2 parameters
4. Compare with single-app staggered results

### Phase 4: Optimization
1. Implement adaptive time stepping in both apps
2. Add sub-cycling capability for friction
3. Optimize transfer efficiency
4. Profile and tune performance

### Phase 5: Extensions
1. 3D implementation (3D bulk + 2D fault SubApp)
2. Multiple fault support (multiple SubApps)
3. Fully dynamic mode with Newmark-Beta
4. Discrete Green's function mode

---

## 9. Input File Templates

### 9.1 Complete MainApp (`multiapp_bp2_main.i`)

```
# SEAS BP2 - Main Elasticity App
# Uses MultiApp for friction

mu = 32.04e9
rho = 2670.0
cs = 3464.0
Vp = 1e-9

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 125
    ny = 50
    xmin = -100000
    xmax = 0
    ymin = 0
    ymax = 40000
    boundary_name_prefix = left
  []
  [right_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 125
    ny = 50
    xmin = 0
    xmax = 100000
    ymin = 0
    ymax = 40000
    boundary_name_prefix = right
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'left_block right_block'
    stitch_boundaries_pairs = 'left_right right_left'
  []
  [fault]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = 0
    paired_block = 1
    new_boundary = fault
  []
[]

[Variables]
  [w]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxVariables]
  [slip_from_friction]
    order = CONSTANT
    family = MONOMIAL
  []
  [traction_to_friction]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = w
  []
[]

[DGKernels]
  [dg_elasticity]
    type = DGElasticityAntiplane
    variable = w
    epsilon = 1.0
    sigma = 6.0
  []
[]

[InterfaceKernels]
  [fault_slip]
    type = DGPrescribedSlipInterfaceKernel
    variable = w
    neighbor_var = w
    boundary = fault
    slip_prescribed = slip_from_friction
    penalty = 1e12
  []
[]

[AuxKernels]
  [compute_traction]
    type = DGElasticTractionAux
    variable = traction_to_friction
    displacement = w
    shear_modulus = ${mu}
    boundary = fault
    execute_on = 'TIMESTEP_END'
  []
[]

[BCs]
  [left_far]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = '-0.5*${Vp}*t'
  []
  [right_far]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = '0.5*${Vp}*t'
  []
  [free_surface_left]
    type = DGElasticityNeumannBC
    variable = w
    boundary = left_bottom
    traction = 0
  []
  [free_surface_right]
    type = DGElasticityNeumannBC
    variable = w
    boundary = right_bottom
    traction = 0
  []
[]

[MultiApps]
  [friction_app]
    type = TransientMultiApp
    input_files = 'multiapp_bp2_friction.i'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Transfers]
  [send_traction]
    type = MultiAppNearestNodeTransfer
    to_multi_app = friction_app
    source_variable = traction_to_friction
    variable = traction_received
    source_boundary = fault
  []
  [receive_slip]
    type = MultiAppNearestNodeTransfer
    from_multi_app = friction_app
    source_variable = slip
    variable = slip_from_friction
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1e6
  num_steps = 100
[]

[Outputs]
  exodus = true
  csv = true
[]
```

### 9.2 Complete SubApp (`multiapp_bp2_friction.i`)

```
# SEAS BP2 - Friction SubApp
# 1D mesh along fault (z-direction)

a = 0.025
b = 0.015
Dc = 0.004
f0 = 0.6
V0 = 1e-6
sigma_n = 50e6
eta = 4.634e6
tau_pre = 26.546e6
Vinit = 1e-9
theta0 = ${fparse Dc/Vinit}

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 50
  xmin = 0
  xmax = 40000  # Fault depth Wf = 40 km
[]

[Variables]
  [state_variable]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${theta0}
  []
[]

[AuxVariables]
  [traction_received]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${tau_pre}
  []
  [slip_rate]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${Vinit}
  []
  [slip]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[Kernels]
  [state_time]
    type = TimeDerivative
    variable = state_variable
  []
  [state_evolution]
    type = StateEvolutionKernel
    variable = state_variable
    slip_rate = slip_rate
    Dc = ${Dc}
  []
[]

[AuxKernels]
  [solve_V]
    type = SlipRateAux
    variable = slip_rate
    traction = traction_received
    state_variable = state_variable
    a = ${a}
    b = ${b}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    tau_pre = ${tau_pre}
    execute_on = 'TIMESTEP_BEGIN'
  []
  [integrate_slip]
    type = SlipIntegrationAux
    variable = slip
    slip_rate = slip_rate
    execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1e6
[]

[Outputs]
  exodus = true
[]
```

---

## 10. Key Differences from Current Implementation

| Current Staggered | MultiApp Proposed |
|-------------------|-------------------|
| AuxKernels on boundary of bulk mesh | Separate 1D SubApp mesh |
| `execute_on = TIMESTEP_BEGIN/END` | Transfer + SubApp solve |
| All variables in one system | Variables split across apps |
| Single `Transient` solve | Nested `TransientMultiApp` |
| `InterfaceMaterial` for traction | `Transfer` for data exchange |

---

## 11. Migration Path

To convert existing single-app code to MultiApp:

1. **Extract friction logic** from AuxKernels to SubApp Kernels
2. **Create 1D fault mesh** matching boundary nodes
3. **Replace InterfaceMaterial** with Transfer
4. **Move friction parameters** to SubApp input
5. **Test incrementally** - first with constant slip, then full coupling

---

## 12. References

1. Tandem Paper: Uphoff et al., "A discontinuous Galerkin method for sequences of earthquakes and aseismic slip on multiple faults using unstructured curvilinear grids"
2. Tandem Code: https://github.com/TEAR-ERC/tandem
3. MOOSE MultiApp Documentation: https://mooseframework.inl.gov/syntax/MultiApps/
4. MOOSE Transfers Documentation: https://mooseframework.inl.gov/syntax/Transfers/
5. SCEC SEAS Benchmark: https://strike.scec.org/cvws/seas/

---

*Document Version: 2.0*
*Updated: January 2026*
*Architecture: MultiApp SEAS Implementation*
