# Legacy regression: exercises the pre-existing darcy_poiseuille_permeability_model
# branch (unchanged by this refactor).
#
# Uniform damage d = 0.5, zero strain (undeformed solid).
# Legacy formula: w = d * wc, k_f = w^2/12, k_eff = k0 + d^n * (k_f - k0)
# then rotated by R (which is a no-op for an isotropic tensor).
#
# With wc = 1e-6, d = 0.5, perm_exponent = 10, intrinsic = 5e-19:
#   w = 5e-7,  k_f = 2.08333e-14,  d^n = 0.5^10 = 9.7656e-4,
#   k_eff = 5e-19 + 9.7656e-4 * (2.08333e-14 - 5e-19) ~= 2.034e-17.
# All three diagonals equal, off-diagonals zero.

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
    function = '0'
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
    function = '0.5'
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
    # Legacy path: enable via the existing boolean flag; the enum stays at "none"
    # so the constructor falls back to the bool.
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = 5e-19
    wc = 1e-6
    perm_exponent = 10
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
