# Method of Manufactured Solutions for 2D Dynamic Simulation with Continuous Domain
# Using exponential function for temporal evolution
[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 101
    ny = 101
    xmin = -2000
    xmax = 2000
    ymin = -2000
    ymax = 2000
    elem_type = QUAD4
  [../]
  
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
  [../]
  [./accel_x]
  [../]
  [./vel_y]
  [../]
  [./accel_y]
  [../]
  [./mms_exact_disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mms_exact_disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./error_field_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./error_field_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]
[Functions]
  # ================ DISPLACEMENT FUNCTIONS ================
  [./mms_displacement_x]
    type = ParsedFunction
    expression = 'delta*exp(-(x^2 + y^2)/(2*width^2) + (t - tbar)/tw)/2'
    symbol_names = 'delta tbar tw width'
    symbol_values = '1.0 0 0.1 100.0'
  [../]

  [./mms_displacement_y]
    type = ParsedFunction
    expression = 'delta*exp(-(x^2 + y^2)/(2*width^2) + (t - tbar)/tw)/2'
    symbol_names = 'delta tbar tw width'
    symbol_values = '1.0 0 0.1 100.0'
  [../]

  # ================ VELOCITY FUNCTIONS ================
  [./initial_velocity_x]
    type = ParsedFunction
    expression = 'delta*exp(-(x^2 + y^2)/(2*width^2) + (t - tbar)/tw)/(2*tw)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '1.0 0 0.1 100.0'
  [../]

  [./initial_velocity_y]
    type = ParsedFunction
    expression = 'delta*exp(-(x^2 + y^2)/(2*width^2) + (t - tbar)/tw)/(2*tw)'
    symbol_names = 'delta tbar tw width'
    symbol_values = '1.0 0 0.1 100.0'
  [../]

  # ================ SOURCE TERMS ================
  [./source_x]
    type = ParsedFunction
    expression = 'delta*(-1.0*mu*tw^2*(-width^2 + y*(x + y)) + rho*width^4 + tw^2*width^2*(lambda + 2*mu) - tw^2*x*(lambda*y + x*(lambda + 2*mu)))*exp((-tw*(x^2 + y^2) + 2*width^2*(t - tbar))/(2*tw*width^2))/(2*tw^2*width^4)'
    symbol_names = 'delta tbar tw width lambda mu rho'
    symbol_values = '1.0 0 0.1 100.0 20e9 30e9 2850'
  [../]

  [./source_y]
    type = ParsedFunction
    expression = 'delta*(-1.0*mu*tw^2*(-width^2 + x*(x + y)) + rho*width^4 + tw^2*width^2*(lambda + 2*mu) - tw^2*y*(lambda*x + y*(lambda + 2*mu)))*exp((-tw*(x^2 + y^2) + 2*width^2*(t - tbar))/(2*tw*width^2))/(2*tw^2*width^4)'
    symbol_names = 'delta tbar tw width lambda mu rho'
    symbol_values = '1.0 0 0.1 100.0 20e9 30e9 2850'
  [../]
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = false
        planar_formulation = PLANE_STRAIN
        generate_output = 'stress_xx stress_yy stress_xy'
      [../]
    [../]
  [../]
[]

[AuxKernels]
  [./velocity_x]
    type = CompVarRate
    variable = vel_x
    coupled = disp_x
  [../]
  
  [./velocity_y]
    type = CompVarRate
    variable = vel_y
    coupled = disp_y
  [../]
  
  [./mms_solution_x]
    type = FunctionAux
    function = mms_displacement_x
    variable = mms_exact_disp_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./mms_solution_y]
    type = FunctionAux
    function = mms_displacement_y
    variable = mms_exact_disp_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./error_field_x_aux]
    type = ParsedAux
    variable = error_field_x
    coupled_variables = 'disp_x mms_exact_disp_x'
    expression = 'abs(disp_x - mms_exact_disp_x)'
    execute_on = 'TIMESTEP_END'
  [../]
  
  [./error_field_y_aux]
    type = ParsedAux
    variable = error_field_y
    coupled_variables = 'disp_y mms_exact_disp_y'
    expression = 'abs(disp_y - mms_exact_disp_y)'
    execute_on = 'TIMESTEP_END'
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
  
  [./source_force_x]
    type = BodyForce
    variable = disp_x
    function = source_x
  [../]
  
  [./source_force_y]
    type = BodyForce
    variable = disp_y
    function = source_y
  [../]
[]


[ICs]
  [./disp_x_ic]
    type = FunctionIC
    variable = disp_x
    function = mms_displacement_x
  [../]
  
  [./disp_y_ic]
    type = FunctionIC
    variable = disp_y
    function = mms_displacement_y
  [../]
  
  [./vel_x_ic]
    type = FunctionIC
    variable = vel_x
    function = initial_velocity_x
  [../]
  
  [./vel_y_ic]
    type = FunctionIC
    variable = vel_y
    function = initial_velocity_y
  [../]
[]

[Materials]
  [./elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 20e9
    shear_modulus = 30e9
    use_displaced_mesh = false
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2850
  [../]
[]

[Postprocessors]
  [./l2_error_x]
    type = ElementL2Error
    function = mms_displacement_x
    variable = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./l2_error_y]
    type = ElementL2Error
    function = mms_displacement_y
    variable = disp_y
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./l2_error_total_x]
    type = ElementL2Norm
    variable = error_field_x
    execute_on = 'TIMESTEP_END'
  [../]
  
  [./l2_error_total_y]
    type = ElementL2Norm
    variable = error_field_y
    execute_on = 'TIMESTEP_END'
  [../]
  
  [./max_error_x]
    type = ElementExtremeValue
    variable = error_field_x
    value_type = max
    execute_on = 'TIMESTEP_END'
  [../]
  
  [./max_error_y]
    type = ElementExtremeValue
    variable = error_field_y
    value_type = max
    execute_on = 'TIMESTEP_END'
  [../]
  
  [./max_vel_x]
    type = ElementExtremeValue
    variable = vel_x
  [../]
  
  [./max_vel_y]
    type = ElementExtremeValue
    variable = vel_y
  [../]
[]

[BCs]
  [./x_bc]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right top bottom'
    function = mms_displacement_x
  [../]
  
  [./y_bc]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'left right top bottom'
    function = mms_displacement_y
  [../]
[]

[Preconditioning]
    [./smp]
        type = SMP
        full = true

    [../]
[]


[Executioner]
  type = Transient
  dt = 0.0001
  end_time = 0.3
  
  [./TimeIntegrator]
    type = CentralDifference
  [../]
[]

[Outputs]
  exodus = true
  time_step_interval = 10
  checkpoint = true
  [./csv]
    type = CSV
    time_data = true
  [../]
[]