[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 200
        ny = 200
        xmin = -4
        xmax = 4
        ymin = -4
        ymax = 4
    [../]
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
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

    [lambda_y]
        # Lagrange multiplier for normal traction
        family = LAGRANGE
        order = FIRST
        initial_condition = 1e-20
    []
    [lambda_x]
        # Lagrange multiplier for tangential traction
        family = LAGRANGE
        order = FIRST
        initial_condition = 1e-20
    []


[]

[AuxVariables]
    # Velocity calculations
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    # Variables to receive slip from subapp
    [slip_x]
        family = MONOMIAL
        order = FIRST
    []
    [slip_xx]
        family = LAGRANGE
        order = FIRST
    []
    [slip_y]
        family = MONOMIAL
        order = FIRST
    []
    [sliprate]
        family = MONOMIAL
        order = FIRST
    []
    [min_elem_size]
         order = CONSTANT
        family = MONOMIAL
    []

[]

[AuxKernels]
    # Calculate velocity
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [min_elem_size_calculator]
        type = ElementLengthAux
        variable = min_elem_size
        method = min
        execute_on = initial
    []
    [base_nodal_proj_lagrange]
        type = ProjectionAux
        variable = slip_xx
        v = slip_x
    []
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_x disp_y'
            planar_formulation = PLANE_STRAIN
            generate_output = 'stress_xx stress_yy stress_xy'
        [../]
        [../]
    [../]
[]

[Kernels]

    # Lagrange multiplier time derivatives
    [lambda_x_time]
        type = CoefTimeDerivative
        variable = lambda_x
        Coefficient = 1e-5 # Small coefficient
    []
    [lambda_y_time]
        type = CoefTimeDerivative
        variable = lambda_y
        Coefficient = 1e-5  # Small coefficient
    []
[]

[Problem]
    kernel_coverage_check = false
[]


[InterfaceKernels]
     [y_slip]
        type = FaultSlipInterfaceKernel
        variable = lambda_y
        neighbor_var = lambda_y
        slip = slip_y
        coupled_disp = disp_y
        boundary = 'Block0_Block1'
    []
    
    [x_slip]
        type = FaultSlipInterfaceKernel
        variable = lambda_x
        neighbor_var = lambda_x
        slip = slip_xx
        coupled_disp = disp_x
        boundary = 'Block0_Block1'
    []
    [lambda_X]
        type = FaultShearTractionLagrangeMultiplier
        variable = disp_x
        neighbor_var = disp_x
        lambda_x = lambda_x
        boundary = 'Block0_Block1'
    []
    
    [lambda_Y]
        type = FaultNormalTractionLagrangeMultiplier
        variable = disp_y
        neighbor_var = disp_y
        lambda_y = lambda_y
        boundary = 'Block0_Block1'
    []
  
    
[]

[Materials]
    # Elasticity and mechanics materials
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 20e9
        shear_modulus = 30e9 
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []

[]


[BCs]
  [roller_xmin]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'left right'
  []
  [roller_ymin]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'bottom top'
  []
  
  
[]



[Postprocessors]
  [max_sliprate_value]
    type = ElementExtremeValue
    variable = sliprate  
    execute_on = 'TIMESTEP_END'
  []
  [min_elem_size_pp]
    type = ElementExtremeValue
    variable = min_elem_size
    value_type = min
    execute_on = 'INITIAL'
  []
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options = '-ksp_view_mat -ksp_view_preconditioner -ksp_monitor_singular_value'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
    petsc_options_value = 'gmres      lu       mumps                        200'
  []
[]

[Executioner]
    type = Transient
 #   solve_type = Newton
    nl_rel_tol = 1e-6  # Relative tolerance
 # nl_abs_tol = 1e-8  # Absolute tolerance
    
    [TimeStepper]
        type = AdaptiveTimeStepCalculator
        cp = 5500                # P-wave velocity
        a_prem = 0.1             # Proportionality index for stable time step
        shear_modulus = 30e9     # Shear modulus
        permeability = 1  # Permeability value (k_py)
        biot_modulus = 1   # Inverse of Biot modulus (c_o)
        a_o = 0.015              # Direct effect parameter
        b_o = 0.02               # State evolution parameter
        L = 1e-5                # Characteristic slip distance
        normal_stress = 50e6    # Normal stress on the fault
        zeta_max = 0.5           # Maximum allowed value for zeta
        dx_min = min_elem_size_pp              # Minimum element size
        max_slip_rate = max_sliprate_value  # Reference to the postprocessor that computes max slip rate
    []
    end_time = 8.0
[]

[Outputs]
  exodus = true
[]




[MultiApps]
    # SubApp to solve rate-and-state friction
    [./rate_state_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'Eric_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    # Send Lagrange multipliers and pressure to subapp
    [send_tractions]
        type = MultiAppCopyTransfer
        to_multi_app = rate_state_app
        source_variable = 'lambda_x lambda_y'
        variable = 'traction_tangential traction_normal'
        execute_on = 'TIMESTEP_BEGIN'
    []
    
    # Get slip from subapp
    [get_slip]
        type = MultiAppCopyTransfer
        from_multi_app = rate_state_app
        source_variable = 'slip_x slip_y sliprate_sub'
        variable = 'slip_x slip_y sliprate'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]