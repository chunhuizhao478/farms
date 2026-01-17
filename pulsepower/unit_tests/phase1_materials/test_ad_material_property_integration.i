# =============================================================================
# Unit Test: AD Material Property Integration with Postprocessors
# =============================================================================
# This test verifies the correct approach for integrating AD material properties
# using AuxVariable + ADMaterialRealAux + ElementIntegralVariablePostprocessor
#
# Note: Using ElementIntegralMaterialProperty directly with AD properties will fail
# with error: "The requested non-AD material property is already retrieved or
# declared as an AD property"
# =============================================================================

E = 50e9
nu = 0.3
K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
  xmin = 0
  xmax = 0.01
  ymin = 0
  ymax = 0.01
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  # AuxVariable to store AD material property for integration
  [psie_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [psie_active_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [damage_ic]
    type = ConstantIC
    variable = d
    value = 0.0
  []
[]

[Kernels]
  [dispkernel_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [dispkernel_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
[]

[AuxKernels]
  # Extract AD material property using ADMaterialRealAux
  [get_psie]
    type = ADMaterialRealAux
    variable = psie_aux
    property = psie
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_psie_active]
    type = ADMaterialRealAux
    variable = psie_active_aux
    property = psie_active
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  [left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [left_y]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  []
  [right_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 1e-6  # Small displacement
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  # AD Phase-field elasticity material (declares AD properties psie, psie_active, etc.)
  [elasticity]
    type = ADSmallDeformationIsotropicElasticityPF
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_inactive = psie_inactive
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    decomposition = VOLDEV
    model_type = AT1
    eta = 1e-6
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1
  num_steps = 1
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
[]

# Correct approach: Use ElementIntegralVariablePostprocessor with AuxVariable
[Postprocessors]
  [total_strain_energy]
    type = ElementIntegralVariablePostprocessor
    variable = psie_aux
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_active_strain_energy]
    type = ElementIntegralVariablePostprocessor
    variable = psie_active_aux
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  csv = true
[]
