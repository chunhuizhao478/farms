[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../meshfile/cylinder_w_hole.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '-0.00177864  0.000445793  0'
    new_boundary = corner_ptr
    input = msh
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
      youngs_modulus = 67e9
      poissons_ratio = 0.32
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

[Functions]
  [cylinder_outer_pressure_x]
    type = CylinderPressureXdirAD
    value = 17.2e6
  []
  [cylinder_outer_pressure_y]
    type = CylinderPressureYdirAD
    value = 17.2e6
  []
[]

[BCs]
  [Pressure_outer_x]
    type = Pressure
    boundary = 'outer'
    variable = disp_x
    function = cylinder_outer_pressure_x
  []
  [Pressure_outer_y]
    type = Pressure
    boundary = 'outer'
    variable = disp_y
    function = cylinder_outer_pressure_y
  []
  [./fix_cptr_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr
    value = 0
  []
  [./fix_cptr_y]
    type = DirichletBC
    variable = disp_y
    boundary = corner_ptr
    value = 0
  []
  [./fix_cptr_z]
    type = DirichletBC
    variable = disp_z
    boundary = corner_ptr
    value = 0
  []
[]