[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 40  # Number of elements in x-direction (2 m / 0.05 m per element)
    ny = 20  # Number of elements in y-direction (1 m / 0.05 m per element)
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 1
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
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    # [accel_x]
    #     type = NewmarkAccelAux
    #     variable = accel_x
    #     displacement = disp_x
    #     velocity = vel_x
    #     beta = 0.25
    #     execute_on = 'TIMESTEP_END'
    # []
    # [vel_x]
    #     type = NewmarkVelAux
    #     variable = vel_x
    #     acceleration = accel_x
    #     gamma = 0.5
    #     execute_on = 'TIMESTEP_END'
    # []
    # [accel_y]
    #     type = NewmarkAccelAux
    #     variable = accel_y
    #     displacement = disp_y
    #     velocity = vel_y
    #     beta = 0.25
    #     execute_on = 'TIMESTEP_END'
    # []
    # [vel_y]
    #     type = NewmarkVelAux
    #     variable = vel_y
    #     acceleration = accel_y
    #     gamma = 0.5
    #     execute_on = 'TIMESTEP_END'
    # []
    [vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
[]

[Kernels]
    [dispkernel_x]
        type = TotalLagrangianStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
    [dispkernel_y]
        type = TotalLagrangianStressDivergence
        variable = disp_y
        component = 1
        large_kinematics = true
    []
    [./inertia_x]
        type = InertialForce
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
    [./inertia_y]
        type = InertialForce
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
    [damping_x]
        type = StiffPropDampingImplicit
        variable = disp_x
        component = 0
        zeta = 0.1
    []
    [damping_y]
        type = StiffPropDampingImplicit
        variable = disp_y
        component = 1
        zeta = 0.1
    [] 
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
    []
[]  

[Functions]
    # [func_pulse]
    #     type = ParsedFunction
    #     # Gaussian pulse: A * exp(-(t-t0)^2/(2*sigma_t^2)) * exp(-(x-x0)^2/(2*sigma_x^2))
    #     expression = 'if(t < t_end, A * exp(-(t-t0)^2/(2*sigma_t^2)) * exp(-(x-x0)^2/(2*sigma_x^2)), 0)'
    #     symbol_names = 'A t0 sigma_t sigma_x x0 t_end'
    #     symbol_values = '1e3 0 1e-6 100 0 1e-5'  # Amplitude, peak time (t0=0), time width, space width, center position, end time
    # []
    [./func_pulse]
        type = ParsedFunction
        # Exponentially decaying pulse: A * exp(-t/tau)
        expression = 'if(t < t_end, A * exp(-t/tau), 0)'
        symbol_names = 'A tau t_end'
        symbol_values = '1e3 1e-6 1e-3'  # Amplitude, decay constant, end time
    [../]
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
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 1e100
    num_steps = 1000
    l_max_its = 15
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 5
    nl_abs_tol = 1e-8
    # petsc_options_iname = '-ksp_type -pc_type'
    # petsc_options_value = 'gmres     hypre'
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    # automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    dt = 1e-6
    # [./TimeIntegrator]
    #     type = NewmarkBeta
    #     beta = 0.25
    #     gamma = 0.5
    # [../]
    [./TimeIntegrator]
        type = FarmsNewmarkBeta
        beta = 0.25
        gamma = 0.5
        factor = 0.9
        threshold = 1e-5
    [../]
[]

[BCs]
    [./pulse_load]
        type = FunctionNeumannBC
        variable = disp_x  # or disp_y depending on direction
        boundary = 'right'   # boundary where load is applied
        function = func_pulse
    [../]
    [fix_left_x]
        type = DirichletBC
        variable = disp_x
        boundary = left
        value = 0
    []
    [fix_left_y]
        type = DirichletBC
        variable = disp_y
        boundary = left
        value = 0
    []
[]

[Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 1
    [../]
[]