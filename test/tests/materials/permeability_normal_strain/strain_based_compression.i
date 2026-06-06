# Unit test: strain-based crack normal under in-plane COMPRESSION (Liu 2024 eqs.
# 29-30) -- regression test for REVIEW.md R-001.
#
# Purpose: a damaged element (d >= threshold) under in-plane compression must NOT
# receive a spurious fracture permeability. In 2D plane strain eps_zz = 0, so when
# both in-plane principal strains are negative the MAXIMUM principal strain is the
# out-of-plane eps_zz = 0 and e_1 = e_z. Without the eps_1 > 0 gate the model would
# set n_F = e_z, eps_nn = 0, w_c = h_c (full aperture), projector diag(1,1,0), and
# enhance the in-plane permeability of a CLOSED crack. The gate
# `have_normal = (eps_1 > 0)` routes this point to matrix perm k0*I instead.
#
# Geometry: single QUAD4 of size 1 x 1 m on [0,1]^2.
# Damage:   uniform d = 0.7 (>= threshold 0.5).
# Strain:   disp_x = -1e-3*x, disp_y = -1e-3*y -> eps = diag(-1e-3, -1e-3, 0).
#           Eigenvalues ascending {-1e-3, -1e-3, 0}; max principal strain eps_1 = 0
#           (the out-of-plane direction). eps_1 <= 0 => no tensile opening.
#
# Expected (with the R-001 fix): K = k0*I (no enhancement).
#   K_xx = K_yy = K_zz = k0 = 5e-19,  K_xy = 0.
#
# Without the fix (legacy, pre-R-001): n_F = e_z, w_c = h_c = 1e-3,
#   k_w = (1e-3)^2/12 = 8.333333e-8, projector diag(1,1,0), d^b = 0.49, so
#   K_xx = K_yy = 0.49*k_w = 4.083333e-8 (SPURIOUS), K_zz = k0. The collapse of
#   K_xx, K_yy from 4.08e-8 back to k0 = 5e-19 is what this test pins down.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [d]
    order = FIRST
    family = LAGRANGE
  []
  [eff_perm_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_22]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_01]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [disp_x_aux]
    type = FunctionAux
    variable = disp_x
    function = '-1e-3 * x'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y_aux]
    type = FunctionAux
    variable = disp_y
    function = '-1e-3 * y'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0.7'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [eff_perm_00]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 0
    variable = eff_perm_00
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_11]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 1
    column = 1
    variable = eff_perm_11
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_22]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 2
    column = 2
    variable = eff_perm_22
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_01]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 1
    variable = eff_perm_01
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [u_left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [u_right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 50e9
    poissons_ratio = 0.3
  []
  [strain]
    type = ComputeSmallStrain
  []
  [bulk]
    type = GenericConstantMaterial
    prop_names = 'K G'
    prop_values = '4.1666666667e10 1.9230769231e10'
  []
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    decomposition = SPECTRAL
    model_type = AT1
    eta = 1e-6
    porous_flow_coupling = true
    permeability_model = normal_strain
    intrinsic_permeability = 5e-19
    perm_exponent = 2
    crack_normal_source = principal_strain   # n_F = e_1 (Liu 2024 eqs. 29-30)
    characteristic_length_type = constant
    characteristic_length_value = 1e-3
    permeability_anisotropic = true
    damage_threshold_for_permeability = 0.5
    correction_factor_fc = 1.0
  []
  [stress]
    type = NDComputeSmallDeformationStress
    elasticity_model = elasticity
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
[]
