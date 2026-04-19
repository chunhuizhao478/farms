# Unit test: Heider (2021) normal-strain permeability, 45-deg rotation.
# This is the critical test for the tangential projector in the global frame.
#
# Damage profile:  d(x,y) = 0.5 + 0.2*(x + y)  -> grad(d) = (0.2, 0.2, 0)
#                  -> n_d = (1,1,0)/sqrt(2).
# Strain:          disp_x = 1e-3 * (x + y), disp_y = 1e-3 * (x + y)
#                  -> eps_xx = eps_yy = 1e-3,  eps_xy = 1e-3.
#                  eps_nn = n_d . eps . n_d = 0.5*(eps_xx + 2*eps_xy + eps_yy)
#                          = 0.5*(1e-3 + 2e-3 + 1e-3) = 2e-3.
# h_c = 1e-3, w_c = 1e-3 * 1.002, k_w = w_c^2 / 12.
# Projector I - n_d (x) n_d =
#   [[1/2, -1/2, 0],
#    [-1/2, 1/2, 0],
#    [0,    0,   1]].
# Let alpha = d_qp^b * k_w. At the QP d_qp ~= 0.5 + 0.2*(0.5+0.5) = 0.7, so
# alpha ~= 0.7^2 * k_w.
# Expected:
#   K_xx = K_yy = k0 + alpha/2
#   K_xy = -alpha/2
#   K_zz = k0 + alpha
#   K_xz = K_yz = 0.

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
    function = '1e-3 * (x + y)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y_aux]
    type = FunctionAux
    variable = disp_y
    function = '1e-3 * (x + y)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0.5 + 0.2 * (x + y)'
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
    crack_normal_source = damage_gradient
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

[Postprocessors]
  # Element-averaged effective_perm components for analytical sign-test CSV.
  [K_xx]
    type = ElementAverageValue
    variable = eff_perm_00
    execute_on = 'TIMESTEP_END'
  []
  [K_yy]
    type = ElementAverageValue
    variable = eff_perm_11
    execute_on = 'TIMESTEP_END'
  []
  [K_zz]
    type = ElementAverageValue
    variable = eff_perm_22
    execute_on = 'TIMESTEP_END'
  []
  [K_xy]
    type = ElementAverageValue
    variable = eff_perm_01
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
