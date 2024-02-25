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
      file =  '../../../../meshgenerator/cdbm/borehole_3d/borehole_example.msh'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
  [./stress_xz_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./stress_zz_saved]
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
  [./strain_xz_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./strain_zz_saved]
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
        generate_output = 'strain_xx strain_xy strain_yy strain_xz strain_yz strain_zz'
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
  [get_stress_xz]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 2
    variable = stress_xz_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_yz]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 2
    variable = stress_yz_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_zz]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 2
    j = 2
    variable = stress_zz_saved
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
  [get_strain_xz]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 0
    j = 2
    variable = strain_xz_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_yz]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 1
    j = 2
    variable = strain_yz_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_zz]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 2
    j = 2
    variable = strain_zz_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_xi]
    type = CompXi3D
    strainxx = strain_xx_saved
    strainxy = strain_xy_saved
    strainyy = strain_yy_saved
    strainxz = strain_xz_saved
    strainyz = strain_yz_saved
    strainzz = strain_zz_saved    
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
      youngs_modulus = 48.5e9
      poissons_ratio = 0.22
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
  # nemesis = true 
  interval = 1
[]

[BCs]
  #
  [disp_top]
    type = DirichletBC
    boundary = top
    variable = disp_y
    value = 0.00000585
  []
  [disp_bottom]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = -0.00000585
  []
  #
  [disp_left]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0000069
  []
  [disp_right]
    type = DirichletBC
    boundary = right
    variable = disp_x
    value = -0.0000069
  []
  #
  [disp_front]
    type = DirichletBC
    boundary = front
    variable = disp_z
    value = -0.0000069
  []
  [disp_back]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0000069
  []
[]
