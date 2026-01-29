# DG Elasticity Antiplane Test
# Tests SIPG formulation for 2D antiplane shear elasticity
# Problem: -mu * Laplacian(w) = 0
# With Dirichlet BC: w = 0 at left, w = 1 at right
# Expected: linear solution w = x/L

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

[Variables]
  [w]
    order = FIRST
    family = MONOMIAL
  []
[]

[Kernels]
  # Volume term for Laplacian: integral of mu*grad(w).grad(v)
  [diffusion]
    type = Diffusion
    variable = w
  []
[]

[DGKernels]
  [dg_elasticity]
    type = DGElasticityAntiplane
    variable = w
    epsilon = 1.0  # SIPG
    sigma = 6.0
    shear_modulus = shear_modulus
  []
[]

[Materials]
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '1.0'
  []
[]

[BCs]
  [left]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left
    value = 0
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [right]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right
    value = 1
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Postprocessors]
  # Check solution at midpoint (should be 0.5)
  [w_midpoint]
    type = PointValue
    variable = w
    point = '0.5 0.5 0'
  []
  # L2 error (if manufactured solution available)
  [l2_error]
    type = ElementL2Error
    variable = w
    function = 'x'
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
