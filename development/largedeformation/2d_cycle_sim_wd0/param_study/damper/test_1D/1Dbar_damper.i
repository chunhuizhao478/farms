[Mesh]
    [msh]
        type = GeneratedMeshGenerator
        dim = 1
        nx = 1000
        xmin = 0
        xmax = 10
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '5 0 0; 4.99 0 0; 5.01 0 0; 4.98 0 0; 5.02 0 0;
                 4.97 0 0; 5.03 0 0; 4.96 0 0; 5.04 0 0; 4.95 0 0;
                 5.05 0 0; 4.94 0 0; 5.06 0 0; 4.93 0 0; 5.07 0 0;
                 4.92 0 0; 5.08 0 0; 4.91 0 0; 5.09 0 0; 4.90 0 0;
                 5.10 0 0; 4.89 0 0; 5.11 0 0; 4.88 0 0; 5.12 0 0;
                 4.87 0 0; 5.13 0 0; 4.86 0 0; 5.14 0 0; 4.85 0 0;
                 5.15 0 0; 4.84 0 0; 5.16 0 0; 4.83 0 0; 5.17 0 0;
                 4.82 0 0; 5.18 0 0; 4.81 0 0; 5.19 0 0; 4.80 0 0;
                 5.20 0 0; 4.79 0 0; 5.21 0 0; 4.78 0 0; 5.22 0 0;
                 4.77 0 0; 5.23 0 0; 4.76 0 0; 5.24 0 0; 4.75 0 0;
                 5.25 0 0; 4.74 0 0; 5.26 0 0; 4.73 0 0; 5.27 0 0;
                 4.72 0 0; 5.28 0 0; 4.71 0 0; 5.29 0 0;'
        new_boundary = corner_ptr
        input = msh
    []
    construct_side_list_from_node_list=true
[]

[GlobalParams]
    displacements = 'disp_x'
    use_displaced_mesh = true
[]


[Variables]
    [disp_x]
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
[]

[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
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
    [./inertia_x]
        type = InertialForce
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2400'
    []
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
    []
    [compute_elasticity]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 25e9
        poissons_ratio = 0.2
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
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 1e100
    num_steps = 1000
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 20
    nl_abs_tol = 1e-12
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    dt = 3e-6
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
        inactive_tsteps = 1
    [../]
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */accel_x */inertia_x */absorbing_right_x */applied_dynamic_loading_at_point'
    # disable_objects = '*/vel_x */accel_x */inertia_x'
    start_time = -1e-12
    end_time = 3e-6 # dt used in the simulation
  []
[../]

[Functions]
    [func_bc]
        type = ParsedFunction
        expression = 'if(t < t_load, 1.8e-7 + 1e-8, 1.8e-7)'
        symbol_names = 't_load'
        symbol_values = '0.02'
    []
    # [./point_load]
    #     type = ParsedFunction
    #     expression = 'if(t < t_load, F0, 0)'
    #     symbol_names = 't_load F0'
    #     symbol_values = '6e-6 1e3'
    # [../]
    [./point_load]
        type = ParsedFunction
        expression = 'if(t <= t_load, F0 * sin(pi * t / t_load), 0)'
        # expression = 'if(t <= t_load, F0 * 0.5 * (1 - cos(2 * pi * t / t_load)), 0)'
        symbol_names = 't_load F0'
        symbol_values = '5e-5 1e2'  # Adjust t_load and F0 as needed
    [../]
    # [./point_load]
    #     type = ParsedFunction
    #     expression = 'F0 * exp(-((t - t0)^2) / (2 * sigma^2))'
    #     symbol_names = 'F0 t0 sigma'
    #     symbol_values = '1e2 0.5 0.2'  # Adjust t0 (peak time) and sigma (width)
    # [../]
[]

[Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 1
    [../]
[]

[BCs]
    [applied_dynamic_loading_at_point]
        type = FunctionNeumannBC
        variable = disp_x
        function = point_load
        boundary = corner_ptr
    []
    [fix_left_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = left
    []
    [applied_force_right_x]
        type = NeumannBC
        variable = disp_x
        value = 1e3
        boundary = right
    []
    [absorbing_right_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x'
        velocities = 'vel_x'
        accelerations = 'accel_x'
        component = 0
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 2082.5
        p_wave_speed = 3403.1 
        density = 2400
    []
[]