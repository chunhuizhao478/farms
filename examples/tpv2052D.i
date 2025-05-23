[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
    xmin = -2000
    xmax = 2000
    ymin = -2000
    ymax = 2000
    elem_type = QUAD4
  []
  [./new_block]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'y>0'
    block_id = 1
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
    split_interface = true
    add_interface_on_two_sides = true
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  use_displaced_mesh = false
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
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
  [./nodal_area]
        order = FIRST
        family = LAGRANGE
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
  [./source_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./source_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  # MMS Displacement Functions for Block 0
  [./mms_displacement_x_block0]
    type = ParsedFunction
    expression = '-(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./mms_displacement_y_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]
  
  # MMS Displacement Functions for Block 1
  [./mms_displacement_x_block1]
    type = ParsedFunction
    expression = '(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./mms_displacement_y_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]
  
  # Initial Velocity Functions for Block 0
  [./initial_velocity_x_block0]
    type = ParsedFunction
    expression = '-(delta/2) * (1/tw) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./initial_velocity_y_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]
  
  # Initial Velocity Functions for Block 1
  [./initial_velocity_x_block1]
    type = ParsedFunction
    expression = '(delta/2) * (1/tw) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./initial_velocity_y_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]
  
  # Derivatives for Block 0
  [./dx_mms_disp_x_block0]
    type = ParsedFunction
    expression = '-(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) * (-x/width^2 - c/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./d2x_mms_disp_x_block0]
    type = ParsedFunction
    expression = '-(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) * (-1/width^2 + x*x/width^4 + 2*c*x/(width^2*width) + c*c/width^2)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  # Derivatives for Block 1
  [./dx_mms_disp_x_block1]
    type = ParsedFunction
    expression = '(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) * (-x/width^2 - c/width)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  [./d2x_mms_disp_x_block1]
    type = ParsedFunction
    expression = '(delta/2) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) * (-1/width^2 + x*x/width^4 + 2*c*x/(width^2*width) + c*c/width^2)'
    symbol_names = 'delta tbar tw width c'
    symbol_values = '1.0 0.25 0.1 100.0 3400'
  [../]
  
  # Stress Functions for Block 0
  [./stress_xx_block0]
    type = ParsedFunction
    expression = '(lambda + 2*mu) * dx'
    symbol_names = 'lambda mu dx'
    symbol_values = '20e9 30e9 dx_mms_disp_x_block0'
  [../]
  
  # Stress Functions for Block 1
  [./stress_xx_block1]
    type = ParsedFunction
    expression = '(lambda + 2*mu) * dx'
    symbol_names = 'lambda mu dx'
    symbol_values = '20e9 30e9 dx_mms_disp_x_block1'
  [../]
  
  # Source Terms for Block 0
  [./source_x_block0]
    type = ParsedFunction
    expression = 'rho * (
      -(1/tw) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) 
      + (lambda + 2*mu) * d2x
    )'
    symbol_names = 'rho lambda mu tw width tbar c d2x'
    symbol_values = '2850.0 20e9 30e9 0.1 100.0 0.25 3400 d2x_mms_disp_x_block0'
  [../]
  
  [./source_y_block0]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]
  
  # Source Terms for Block 1
  [./source_x_block1]
    type = ParsedFunction
    expression = 'rho * (
      -(1/tw) * exp((t - tbar)/tw) * exp(-x*x/(2*width*width)) * exp(-c*x/width) 
      + (lambda + 2*mu) * d2x
    )'
    symbol_names = 'rho lambda mu tw width tbar c d2x'
    symbol_values = '2850.0 20e9 30e9 0.1 100.0 0.25 3400 d2x_mms_disp_x_block1'
  [../]
  
  [./source_y_block1]
    type = ParsedFunction
    expression = '0'
    symbol_names = 'delta'
    symbol_values = '1.0'
  [../]

[]

[BCs]
  # Neumann Boundary Conditions at Fault Interface for Stress
  [./fault_lower_stress_xx]
    type = NeumannBC
    variable = disp_x
    boundary = 'Block0_Block1'
    function = stress_xx_block0
  [../]
  
  [./fault_lower_stress_xy]
    type = NeumannBC
    variable = disp_y
    boundary = 'Block0_Block1'
    function = stress_xy_block0
  [../]
  
  [./fault_upper_stress_xx]
    type = NeumannBC
    variable = disp_x
    boundary = 'Block1_Block0'
    function = stress_xx_block1
  [../]
  
  [./fault_upper_stress_xy]
    type = NeumannBC
    variable = disp_y
    boundary = 'Block1_Block0'
    function = stress_xy_block1
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
        block = '0 1'
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
  
  # Block-specific MMS solution aux variables
  [./mms_solution_x_block0]
    type = FunctionAux
    function = mms_displacement_x_block0
    variable = mms_exact_disp_x
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./mms_solution_y_block0]
    type = FunctionAux
    function = mms_displacement_y_block0
    variable = mms_exact_disp_y
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./mms_solution_x_block1]
    type = FunctionAux
    function = mms_displacement_x_block1
    variable = mms_exact_disp_x
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./mms_solution_y_block1]
    type = FunctionAux
    function = mms_displacement_y_block1
    variable = mms_exact_disp_y
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  # Block-specific error calculations
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

  [./source_x_aux_block0]
    type = FunctionAux
    function = source_x_block0
    variable = source_x
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./source_y_aux_block0]
    type = FunctionAux
    function = source_y_block0
    variable = source_y
    block = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
    [./source_x_aux_block1]
    type = FunctionAux
    function = source_x_block1
    variable = source_x
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  
  [./source_y_aux_block1]
    type = FunctionAux
    function = source_y_block1
    variable = source_y
    block = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    block = '0 1'
  [../]
  
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    block = '0 1'
  [../]
  
  # Block-specific source forces
  [./source_force_x_block0]
    type = BodyForce
    variable = disp_x
    function = source_x_block0
    block = 0
  [../]
  
  [./source_force_y_block0]
    type = BodyForce
    variable = disp_y
    function = source_y_block0
    block = 0
  [../]
  
  [./source_force_x_block1]
    type = BodyForce
    variable = disp_x
    function = source_x_block1
    block = 1
  [../]
  
  [./source_force_y_block1]
    type = BodyForce
    variable = disp_y
    function = source_y_block1
    block = 1
  [../]
[]

[ICs]
  # Block-specific initial conditions
  [./disp_x_ic_block0]
    type = FunctionIC
    variable = disp_x
    function = mms_displacement_x_block0
    block = 0
  [../]
  
  [./disp_y_ic_block0]
    type = FunctionIC
    variable = disp_y
    function = mms_displacement_y_block0
    block = 0
  [../]
  
  [./vel_x_ic_block0]
    type = FunctionIC
    variable = vel_x
    function = initial_velocity_x_block0
    block = 0
  [../]
  
  [./vel_y_ic_block0]
    type = FunctionIC
    variable = vel_y
    function = initial_velocity_y_block0
    block = 0
  [../]
  
  [./disp_x_ic_block1]
    type = FunctionIC
    variable = disp_x
    function = mms_displacement_x_block1
    block = 1
  [../]
  
  [./disp_y_ic_block1]
    type = FunctionIC
    variable = disp_y
    function = mms_displacement_y_block1
    block = 1
  [../]
  
  [./vel_x_ic_block1]
    type = FunctionIC
    variable = vel_x
    function = initial_velocity_x_block1
    block = 1
  [../]
  
  [./vel_y_ic_block1]
    type = FunctionIC
    variable = vel_y
    function = initial_velocity_y_block1
    block = 1
  [../]
[]

[Materials]
  [./elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 30e9
    shear_modulus = 30e9
    block = '0 1'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0 1'
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2850
    block = '0 1'
  [../]
[]

[UserObjects]
  [./nodal_area]
    type = NodalArea
    variable = nodal_area
    boundary = Block0_Block1
    execute_on = 'initial TIMESTEP_BEGIN'
  [../]
[]


[Postprocessors]
  # Block-specific L2 error calculations
  [./l2_error_x_block0]
    type = ElementL2Error
    function = mms_displacement_x_block0
    variable = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  [../]
  
  [./l2_error_y_block0]
    type = ElementL2Error
    function = mms_displacement_y_block0
    variable = disp_y
    execute_on = 'INITIAL TIMESTEP_END'
    block = 0
  [../]
  
  [./l2_error_x_block1]
    type = ElementL2Error
    function = mms_displacement_x_block1
    variable = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
    block = 1
  [../]
  
  [./l2_error_y_block1]
    type = ElementL2Error
    function = mms_displacement_y_block1
    variable = disp_y
    execute_on = 'INITIAL TIMESTEP_END'
    block = 1
  [../]
  
  # Total L2 norms of the error fields
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
  
  # Maximum errors
  [./max_error_x_block0]
    type = ElementExtremeValue
    variable = error_field_x
    value_type = max
    execute_on = 'TIMESTEP_END'
    block = 0
  [../]
  
  [./max_error_y_block0]
    type = ElementExtremeValue
    variable = error_field_y
    value_type = max
    execute_on = 'TIMESTEP_END'
    block = 0
  [../]
  
  [./max_error_x_block1]
    type = ElementExtremeValue
    variable = error_field_x
    value_type = max
    execute_on = 'TIMESTEP_END'
    block = 1
  [../]
  
  [./max_error_y_block1]
    type = ElementExtremeValue
    variable = error_field_y
    value_type = max
    execute_on = 'TIMESTEP_END'
    block = 1
  [../]
  
  # Velocity maximums
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
  # External boundaries - use block-specific functions
  [./x_bc_block0]
    type = MatchedValueBC
    variable = disp_x
    boundary = 'left right bottom'
    v = mms_exact_disp_x
  [../]
  
  [./y_bc_block0]
    type = MatchedValueBC
    variable = disp_y
    boundary = 'left right bottom'
    v = mms_exact_disp_y
  [../]
  
  [./x_bc_block1]
    type = MatchedValueBC
    variable = disp_x
    boundary = 'left right top'
     v = mms_exact_disp_x
  [../]
  
  [./y_bc_block1]
    type = MatchedValueBC
    variable = disp_y
    boundary = 'left right top'
     v = mms_exact_disp_y
  [../]


[]

[Executioner]
  type = Transient
  dt = 0.0001
  end_time = 0.3
  
  [./TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  [../]
  petsc_options = '-snes_converged_reason'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = 'lu superlu_dist 100'
  line_search = 'none'
[]

[Outputs]
  exodus = true
  time_step_interval = 20
  [./csv]
    type = CSV
    time_data = true
  [../]
[]