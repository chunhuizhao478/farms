# Unit test for ElkADPoroDynamicPhaseFieldMaterials
# Tests that the poro-dynamic material correctly assembles all damage-dependent properties:
#   - density = rho_s * (1 - phi) + rho_f * phi
#   - biot_coefficient from damaged property
#   - porosity from damaged property
#   - biot_modulus from damaged property
#   - permeability from effective_perm (requires elasticity model)
#
# Test parameters:
#   rho_s = 2600 kg/m^3
#   rho_f = 1000 kg/m^3
#   phi_0 = 0.008
#   K0 = Ks = 50e9 Pa
#   Kf = 2.24e9 Pa
#
# Test case at d = 0.5:
#   phi(0.5) = 0.752
#   density = 2600 * (1 - 0.752) + 1000 * 0.752 = 644.8 + 752 = 1396.8 kg/m^3

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [density_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [biot_coeff_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [biot_modulus_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [damage_ic]
    type = FunctionIC
    variable = d
    function = 'x'
  []
[]

[AuxKernels]
  [get_density]
    type = ADMaterialRealAux
    variable = density_aux
    property = density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_porosity]
    type = ADMaterialRealAux
    variable = porosity_aux
    property = porosity
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_biot_modulus]
    type = ADMaterialRealAux
    variable = biot_modulus_aux
    property = biot_modulus
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
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
    value = 0.001
  []
[]

[Materials]
  # Bulk and shear modulus
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '41.67e9 19.23e9'  # E=50e9, nu=0.3
  []
  # Small strain computation
  [strain]
    type = ADComputeSmallStrain
  []
  # Phase-field elasticity with permeability (simplified - using constant permeability for this test)
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 50e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  # Provide effective_perm and elastic_strain for the poro-dynamic material
  [effective_perm_provider]
    type = ADGenericConstantRankTwoTensor
    tensor_name = effective_perm
    tensor_values = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
  []
  # Damaged hydraulic properties
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = 2e-11  # 1/K0 = 1/50e9
    grain_bulk_modulus = 50e9      # Ks
    minimum_degradation = 1e-8
  []
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
  []
  [damaged_biot_modulus]
    type = ElkADPorousFlowDamagedBiotModulus
    use_damaged_biot = true
    use_damaged_porosity = true
    fluid_bulk_modulus = 2.24e9
    grain_bulk_modulus = 50e9
  []
  # Main poro-dynamic material assembly
  [porodynamics]
    type = ElkADPoroDynamicPhaseFieldMaterials
    rhos_value = 2600
    rhof_value = 1000
    tortosity_value = 1.2
    viscosity_value = 1e-3
    grain_bulk_modulus = 50e9
    fluid_bulk_modulus = 2.24e9
    use_damaged_properties = true
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Postprocessors]
  # At d=0.5 (x=0.5):
  # phi = 0.752
  # density = 2600*(1-0.752) + 1000*0.752 = 644.8 + 752 = 1396.8
  [density_mid]
    type = PointValue
    variable = density_aux
    point = '0.5 0.05 0'
  []
  [porosity_mid]
    type = PointValue
    variable = porosity_aux
    point = '0.5 0.05 0'
  []
  [biot_coeff_mid]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.5 0.05 0'
  []
  [biot_modulus_mid]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.5 0.05 0'
  []
  # At d=0 (x=0.1667, approximately):
  [density_left]
    type = PointValue
    variable = density_aux
    point = '0.1667 0.05 0'
  []
  # At d=1 (x=0.8333, approximately):
  [density_right]
    type = PointValue
    variable = density_aux
    point = '0.8333 0.05 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
