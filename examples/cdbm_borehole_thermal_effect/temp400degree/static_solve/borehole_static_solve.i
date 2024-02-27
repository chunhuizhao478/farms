##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4 [Check!]
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677

# under the stress state, 
# -135MPa 70MPa; 70MPa -120e6
# the principal stress is 198e6 

# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.5
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.4) * 198e6) = 233.6m

# Use mesh size = 25m, resolved by more than 7 elements [Check!]

# Use 1% overstress

# Use dt = 2e-4
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../../meshgenerator/cdbm/borehole/drytest/borehole_wofaults.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0  -0.06  0'
    new_boundary = corner_ptr_bottom
    input = msh
  []
  [./extranodeset2]
      type = ExtraNodesetGenerator
      coord = '0.06  0  0'
      new_boundary = corner_ptr_right
      input = extranodeset1
  []
  [./extranodeset3]
      type = ExtraNodesetGenerator
      coord = '0   0.06  0'
      new_boundary = corner_ptr_top
      input = extranodeset2
  []
  [./extranodeset4]
      type = ExtraNodesetGenerator
      coord = '-0.06  0  0'
      new_boundary = corner_ptr_left
      input = extranodeset3
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./stress_xx_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./stress_yy_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./strain_yy_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xi_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'strain_xx strain_xy strain_yy'
      [../]
    [../]
  [../]
[]

[AuxKernels]
  [get_stress_xx]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 0
    variable = stress_xx_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_xy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 1
    variable = stress_xy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_yy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 1
    variable = stress_yy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  #
  [get_strain_xx]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 0
    j = 0
    variable = strain_xx_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_xy]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 0
    j = 1
    variable = strain_xy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_yy]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 1
    j = 1
    variable = strain_yy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_xi]
    type = CompXi
    strainxx = strain_xx_saved
    strainxy = strain_xy_saved
    strainyy = strain_yy_saved
    variable = xi_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
[]

[Materials]
  [stress_medium]
      type = ComputeLinearElasticStress
  []
  #young's modulus: 6.67GPa
  #poisson's ratio: 0.22
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 35e9
      poissons_ratio = 0.1
  []
[]

[Preconditioning]
  [smp]
    full = true
    type = SMP
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       NONZERO               1e-15'
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  interval = 1
[]

[BCs]
  [Pressure_top]
    type = Pressure
    boundary = top
    variable = disp_y
    factor = 20e6
  []
  [Pressure_bottom]
    type = Pressure
    boundary = bottom
    variable = disp_y
    factor = 20e6
  []
  [Pressure_left]
    type = Pressure
    boundary = left
    variable = disp_x
    factor = 40e6
  []
  [Pressure_right]
    type = Pressure
    boundary = right
    variable = disp_x
    factor = 40e6
  []
  [./fix_cptr1_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr_bottom
    value = 0
  []
  [./fix_cptr2_y]
    type = DirichletBC
    variable = disp_y
    boundary = corner_ptr_right
    value = 0
  []
  [./fix_cptr3_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr_top
    value = 0
  []
  [./fix_cptr4_y]
    type = DirichletBC
    variable = disp_y
    boundary = corner_ptr_left
    value = 0
  []
[]
