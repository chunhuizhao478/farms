# Error test for R-003: regularize_crack_normal = true must be accompanied by an
# explicit crack_normal_regularization. This input is identical to
# regularized_normal_core.i EXCEPT it omits crack_normal_regularization, so the
# constructor must paramError (the 1e-8 default is far below |grad(d)| ~ 1/l and
# would make the regularization a silent no-op). Run as a RunException test
# (no gold file needed).

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
[]

[AuxKernels]
  [disp_x_aux]
    type = FunctionAux
    variable = disp_x
    function = '1e-3 * x'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y_aux]
    type = FunctionAux
    variable = disp_y
    function = '0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0.7'
    execute_on = 'INITIAL TIMESTEP_END'
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
    regularize_crack_normal = true
    # crack_normal_regularization deliberately omitted -> must paramError
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
