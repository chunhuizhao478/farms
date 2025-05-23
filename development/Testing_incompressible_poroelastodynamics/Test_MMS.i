# Wave Propagation with Method of Manufactured Solutions
[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 101
    ny = 101
    xmin = -100
    xmax = 100
    ymin = -100
    ymax = 100
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./vel_x]
    order = FIRST
    family = LAGRANGE
  []
  [./vel_y]
    order = FIRST
    family = LAGRANGE
  []
  [./analytical_x]
    order = FIRST
    family = LAGRANGE
  []
  [./analytical_y]
    order = FIRST
    family = LAGRANGE
  []
  [./error_x]
    order = FIRST
    family = LAGRANGE
  []
  [./error_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [./disp_x_ic]
    type = FunctionIC
    variable = disp_x
    function = mms_x
  [../]
  [./disp_y_ic]
    type = FunctionIC
    variable = disp_y
    function = mms_y
  [../]
[]

[Functions]
  [./mms_x]
    type = AnalyticalSolution
    amplitude = 1.0e-3
    sigma = 30.0
    x0 = 0.0
    y0 = 0.0
    t0 = 0.05
    t_width = 0.02
    omega = 31.4159  # 5 Hz
    component = 0
  [../]
  [./mms_y]
    type = AnalyticalSolution
    amplitude = 1.0e-3
    sigma = 30.0
    x0 = 0.0
    y0 = 0.0
    t0 = 0.05
    t_width = 0.02
    omega = 31.4159  # 5 Hz
    component = 1
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
  [../]
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_y
  [../]
  [./source_x]
    type = GaussianSourceTerm
    variable = disp_x
    amplitude = 1.0e-3
    sigma = 30.0
    x0 = 0.0
    y0 = 0.0
    t0 = 0.05
    t_width = 0.02
    omega = 31.4159  # 5 Hz
    component = 0
  [../]
  [./source_y]
    type = GaussianSourceTerm
    variable = disp_y
    amplitude = 1.0e-3
    sigma = 30.0
    x0 = 0.0
    y0 = 0.0
    t0 = 0.05
    t_width = 0.02
    omega = 31.4159  # 5 Hz
    component = 1
  [../]
[]

[AuxKernels]
  [./vel_x_aux]
    type = CompVarRate
    variable = vel_x
    coupled = disp_x
  []
  [./vel_y_aux]
    type = CompVarRate
    variable = vel_y
    coupled = disp_y
  []
  [./analytical_x_aux]
    type = FunctionAux
    variable = analytical_x
    function = mms_x
    execute_on = 'initial timestep_end'
  []
  [./analytical_y_aux]
    type = FunctionAux
    variable = analytical_y
    function = mms_y
    execute_on = 'initial timestep_end'
  []
  [./error_x_aux]
    type = ParsedAux
    variable = error_x
    coupled_variables = 'disp_x analytical_x'
    expression = 'abs(disp_x - analytical_x)'
    execute_on = 'timestep_end'
  []
  [./error_y_aux]
    type = ParsedAux
    variable = error_y
    coupled_variables = 'disp_y analytical_y'
    expression = 'abs(disp_y - analytical_y)'
    execute_on = 'timestep_end'
  []
[]

[Materials]
  [./elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 1.15e11  # Pa (for steel)
    shear_modulus = 7.69e10  # Pa (for steel)
  []
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  []
  [./stress]
    type = ComputeLinearElasticStress
  []
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = 8000.0  # kg/m³ (for steel)
  []
[]

[BCs]
  [./x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]
  [./y_left]
    type = DirichletBC
    variable = disp_y
    boundary = 'left'
    value = 0
  [../]
  [./x_right]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  [../]
  [./y_right]
    type = DirichletBC
    variable = disp_y
    boundary = 'right'
    value = 0
  [../]
  [./x_top]
    type = DirichletBC
    variable = disp_x
    boundary = 'top'
    value = 0
  [../]
  [./y_top]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0
  [../]
  [./x_bottom]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0
  [../]
  [./y_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
[]

[Postprocessors]
  [./l2_error_x]
    type = ElementL2Error
    variable = disp_x
    function = mms_x
    execute_on = 'timestep_end'
  []
  [./l2_error_y]
    type = ElementL2Error
    variable = disp_y
    function = mms_y
    execute_on = 'timestep_end'
  []
  [./max_error_x]
    type = ElementExtremeValue
    variable = error_x
    value_type = max
    execute_on = 'timestep_end'
  []
  [./max_error_y]
    type = ElementExtremeValue
    variable = error_y
    value_type = max
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = Transient
  dt = 1.0e-4
  end_time = 0.1
  
  [TimeIntegrator]
    type = CentralDifference
  []
  
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  l_max_its = 50
  nl_max_its = 10
[]

[Outputs]
  exodus = true
  interval = 5
  csv = true
[]