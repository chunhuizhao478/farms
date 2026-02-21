# Phase-Field Cohesive Zone Model (PF-CZM) with Finite Deformation

## Motivation

The small-strain PF-CZM formulation fails for this three-point bending problem because the material properties (E = 4.41 MPa, sigma_t = 6.43 MPa) yield a critical strain of approximately 146%. Under the small-strain assumption (epsilon << 1), the model never reaches the strain threshold needed to initiate damage. Switching to a finite deformation (large-strain) framework based on Hencky strain resolves this by properly accounting for geometric nonlinearity at large strains.

---

## 1. Kinematics: Deformation Gradient and Hencky Strain

### Deformation Gradient

The deformation gradient maps material points from the reference to the current configuration:

    F = I + grad(u)

where u is the displacement field. In the finite deformation formulation, the stress equilibrium is solved on the displaced mesh (`use_displaced_mesh = true`).

### Right Cauchy-Green Tensor

    C = F^T F

### Hencky (Logarithmic) Strain

The Hencky strain tensor is defined as:

    epsilon = (1/2) ln(C) = (1/2) ln(F^T F)

This is computed via eigendecomposition of C:

    C = V diag(lambda_i) V^T

    epsilon = V diag((1/2) ln(lambda_i)) V^T

where lambda_i are the eigenvalues of C and V contains the corresponding eigenvectors. The Hencky strain is the natural extension of the infinitesimal strain to finite deformations: for small strains, ln(C)/2 reduces to the symmetric gradient of displacement.

---

## 2. Hyperelastic Strain Energy Density

The Hencky-type hyperelastic energy density for an isotropic material is:

    psi_e = (1/2) K [tr(epsilon)]^2 + G dev(epsilon) : dev(epsilon)

where:
- K = E / [3(1 - 2nu)] is the bulk modulus
- G = E / [2(1 + nu)] is the shear modulus
- tr(epsilon) is the volumetric (trace) part
- dev(epsilon) = epsilon - (1/3) tr(epsilon) I is the deviatoric part

The stress conjugate to the Hencky strain is the **Mandel stress**:

    Sigma = d(psi_e)/d(epsilon) = K tr(epsilon) I + 2G dev(epsilon)

---

## 3. Strain Energy Decomposition

To prevent damage under compression, the strain energy is split into an "active" part (which drives fracture) and an "inactive" part (which remains undegraded). Two decomposition strategies are implemented:

### 3a. Volumetric-Deviatoric (VOLDEV) Decomposition

The VOLDEV split separates the volumetric strain into tensile and compressive contributions:

    epsilon_tr  = tr(epsilon)
    epsilon_tr+ = <epsilon_tr>_+  = max(epsilon_tr, 0)    (tensile volumetric strain)
    epsilon_tr- = epsilon_tr - epsilon_tr+                 (compressive volumetric strain)
    epsilon_dev = dev(epsilon)

**Active energy** (drives fracture, degraded by damage):

    psi_active = (1/2) K (epsilon_tr+)^2 + G epsilon_dev : epsilon_dev

**Inactive energy** (compression, never degraded):

    psi_inactive = (1/2) K (epsilon_tr-)^2

**Total energy with degradation:**

    psi_e = g(d) psi_active + psi_inactive

**Mandel stress decomposition:**

    Sigma_pos = K epsilon_tr+ I + 2G epsilon_dev
    Sigma_neg = K epsilon_tr- I
    Sigma     = g(d) Sigma_pos + Sigma_neg

**Physical interpretation:** Under hydrostatic compression (epsilon_tr < 0), the volumetric energy is not degraded, so the material retains full stiffness in compression. The deviatoric energy is always part of the active energy.

### 3b. Spectral Decomposition

The spectral split decomposes the strain tensor by its principal values, providing a more physically accurate separation for mixed-mode fracture:

**Eigendecomposition of the Hencky strain:**

    epsilon = sum_i epsilon_i (n_i x n_i)

where epsilon_i are principal Hencky strains and n_i are principal directions.

**Positive and negative projections:**

    epsilon+ = sum_i <epsilon_i>_+ (n_i x n_i)      (tensile part)
    epsilon- = sum_i <epsilon_i>_- (n_i x n_i)      (compressive part)

where <x>_+ = max(x, 0) is the Macaulay bracket.

**Lame parameter:**

    lambda = K - 2G/3

**Active energy:**

    psi_active = (1/2) lambda (<tr(epsilon)>_+)^2 + G epsilon+ : epsilon+

**Inactive energy:**

    psi_inactive = (1/2) lambda (<tr(epsilon)>_-)^2 + G epsilon- : epsilon-

**Total energy with degradation:**

    psi_e = g(d) psi_active + psi_inactive

**Mandel stress:**

    Sigma_pos = lambda <tr(epsilon)>_+ I + 2G epsilon+
    Sigma_neg = lambda <tr(epsilon)>_- I + 2G epsilon-
    Sigma     = g(d) Sigma_pos + Sigma_neg

**Physical interpretation:** Each principal direction is treated independently. A principal direction under compression retains full stiffness even when other directions are in tension. This is more accurate for mixed-mode (tension + shear) fracture but requires two eigendecompositions per quadrature point per iteration (one for log(C), one for the spectral split of epsilon).

### 3c. Comparison

| Aspect                    | VOLDEV                          | SPECTRAL                          |
|---------------------------|---------------------------------|-----------------------------------|
| Split criterion           | Sign of volumetric strain       | Sign of each principal strain     |
| Active deviatoric energy  | Always active                   | Only tensile principals active    |
| Computational cost        | Low (trace check only)          | High (eigendecomposition)         |
| Convergence               | Robust                          | May struggle (see Section 7)      |
| Mixed-mode accuracy       | Approximate                     | More physically accurate          |
| Recommended for           | Initial validation              | Final production runs             |

---

## 4. Phase-Field Cohesive Zone Model (PF-CZM)

### Crack Geometric Function

The crack surface density follows the AT2 model:

    alpha(d) = d

with derived constants:
- xi = d(alpha)/dd |_{d=0} = 1
- c_0 = 4 integral_0^1 sqrt(alpha(s)) ds = 8/3

### Rational Degradation Function (PF-CZM specific)

The PF-CZM uses a rational degradation function that embeds the material strength:

    g(d) = [(1-d)^p / ((1-d)^p + a_1 d (1 + a_2 d + a_2 a_3 d^2))] (1 - eta) + eta

where:
- p = 2 (exponent)
- a_1 = Gc / (psi_c xi / (c_0 l)) = Gc c_0 l / (psi_c xi) (couples fracture toughness, critical energy, and length scale)
- a_2 = -0.5
- a_3 = 0
- eta = 1e-6 (residual stiffness to prevent singular stiffness matrix)

This degradation function is specifically designed so that:
1. g(0) = 1 (intact material)
2. g(1) = eta (fully damaged)
3. The stress-strain curve exhibits a peak stress equal to sigma_t before softening
4. The total dissipated energy during complete fracture equals Gc

### Critical Energy Density

    psi_c = sigma_t^2 / (2E)

This is the elastic energy density at the tensile strength under uniaxial stress.

### Fracture Free Energy

The total free energy density in the fracture sub-problem:

    psi = alpha(d) Gc / (c_0 l) + g(d) psi_active

The phase-field evolution is obtained by minimizing psi with respect to d, subject to the irreversibility constraint (d can only increase).

---

## 5. Stress Transformation

The Mandel stress Sigma is work-conjugate to the Hencky strain. To obtain the First Piola-Kirchhoff stress P (used in the momentum equation on the reference configuration):

    P = F Sigma / det(F)

The Cauchy stress (true stress in the current configuration) is:

    sigma = (1/J) F Sigma F^T

where J = det(F).

---

## 6. Multi-App Staggered Solution Strategy

The coupled problem is solved using a staggered (alternate minimization) scheme with two sub-applications:

### Elasticity App (this file)
1. **Input:** phase-field d (from fracture app, previous iteration)
2. **Compute:** deformation gradient F from displacements
3. **Compute:** Hencky strain epsilon = (1/2) ln(F^T F)
4. **Compute:** decomposed stress and psi_active
5. **Solve:** momentum balance div(P) = 0 on displaced mesh
6. **Output:** psi_active to fracture app

### Fracture App (fracture.i, unchanged)
1. **Input:** psi_active from elasticity app
2. **Compute:** total free energy psi = alpha(d) Gc/(c_0 l) + g(d) psi_active
3. **Solve:** phase-field evolution with irreversibility bound
4. **Output:** d to elasticity app

The two apps iterate within each time step (fixed-point iteration, up to 20 iterations) until convergence (relative tolerance 1e-6).

---

## 7. Block Structure

The mesh is divided into two blocks:

- **Block 1 (damage zone):** The main body where fracture can occur. Uses the full PF-CZM material stack with degradation function g(d) and strain energy decomposition.

- **Block 2 (elastic zone):** Regions around the supports and loading points where stress concentrations are non-physical artifacts. Uses a `NoDegradation` function (g = 1, no damage regardless of d) with Hencky elasticity but no decomposition. This prevents spurious damage at the boundary condition application points.

---

## 8. Convergence Considerations

### VOLDEV Version
- More robust convergence due to smooth Jacobian
- The Macaulay bracket on tr(epsilon) introduces only one non-smoothness point
- Recommended for initial validation of the finite deformation setup

### SPECTRAL Version
- Two eigendecompositions per QP per iteration: one for ln(C), one for spectral split
- AD derivative singularities when eigenvalues of epsilon coincide (repeated eigenvalues)
- Non-smooth Jacobian at zero Hencky strain eigenvalues (Macaulay bracket applied to each principal strain)
- May require smaller initial time steps and tighter solver tolerances
- If convergence issues arise, validate with VOLDEV first

### Adaptive Time Stepping
Both versions use `IterationAdaptiveDT` with:
- Initial dt = 0.01
- Target 6 nonlinear iterations per step
- Growth factor 1.2, cutback factor 0.5
- Time step limited by `dt_limit_pp` (caps at 0.001 after t = 16.5, when damage is expected)

### Volumetric Locking Correction
Both versions enable `volumetric_locking_correction = true` in GlobalParams, which applies B-bar correction to prevent volumetric locking in low-order (linear) hexahedral/tetrahedral elements, particularly important for nearly incompressible behavior (nu = 0.25).

---

## 9. Summary of Parameters

| Parameter  | Value       | Description                          |
|------------|-------------|--------------------------------------|
| E          | 4.41 MPa    | Young's modulus                      |
| nu         | 0.25        | Poisson's ratio                      |
| K          | 2.94 MPa    | Bulk modulus                         |
| G          | 1.764 MPa   | Shear modulus                        |
| sigma_t    | 6.43 MPa    | Tensile strength                     |
| Gc         | 3 J/m^2     | Fracture toughness                   |
| psi_c      | 4688 J/m^3  | Critical energy density              |
| l          | 5e-5 m      | Regularization length                |
| Loading    | -1.667e-6 t | Displacement-controlled (m per unit time) |

---

## References

- Tian, F., Tang, X., Xu, T., & Li, L. (2023). An adaptive edge-based smoothed finite element method (ES-FEM) for phase-field modeling of fractures at large deformations. *Computer Methods in Applied Mechanics and Engineering*.
- Wu, J.-Y. (2017). A unified phase-field theory for the mechanics of damage and quasi-brittle failure. *Journal of the Mechanics and Physics of Solids*, 103, 72-99.
- Miehe, C., Welschinger, F., & Hofacker, M. (2010). Thermodynamically consistent phase-field models of fracture. *International Journal for Numerical Methods in Engineering*, 83(10), 1273-1311.
- Raccoon: A phase-field fracture framework (https://github.com/hugary1995/raccoon)
