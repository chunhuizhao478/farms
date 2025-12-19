#implicit continuum damage-breakage model with adaptive dynamic/quasi-dynamic switching
#Switches based on deviatoric strain rate threshold (1e-5 1/s)
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../../mesh/mesh_large.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '-480000 -480000 0;
                 480000 -480000 0'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y'
[]

[GlobalParams]
    displacements = 'disp_x disp_y'

    ##----continuum damage breakage model----##
    lambda_o = 32.04e9
    shear_modulus_o = 32.04e9
    xi_0 = -0.8
    xi_d = -0.9
    C_g = 1e-16
    m1 = 10
    m2 = 1
    chi = 0.8
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
    # Velocities (used in both modes)
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [vel_mag]
        order = FIRST
        family = LAGRANGE
    []
    # Accelerations (only for dynamic mode)
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    # Damage/breakage variables
    [alpha_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    [B_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    # Strain and stress outputs
    [I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [deviatroic_strain_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [structural_stress_coefficient_aux]
        order = FIRST
        family = MONOMIAL
    []
    [gradx_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    [grady_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    [cg_aux]
        order = FIRST
        family = LAGRANGE
    []
    [nonlocal_xi]
        order = FIRST
        family = MONOMIAL
    []
    [pk2_stress_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [green_lagrange_elastic_strain_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [plastic_strain_01]
        order = CONSTANT
        family = MONOMIAL
    []
    [total_lagrange_strain_01]
        order = CONSTANT
        family = MONOMIAL
    []
    ##
    [eqstrain_nonlocal_initial_aux]
        order = FIRST
        family = MONOMIAL
    []
[]

[AuxKernels]
    # Velocity computation (quasi-dynamic mode)
    [vel_x_aux_quasi]
        type = FarmsVelocityAux
        variable = vel_x
        displacement = disp_x
        execute_on = 'TIMESTEP_END'
        enable = true  # Start in quasi-dynamic mode
    []
    [vel_y_aux_quasi]
        type = FarmsVelocityAux
        variable = vel_y
        displacement = disp_y
        execute_on = 'TIMESTEP_END'
        enable = true
    []
    # Velocity computation (dynamic mode using Newmark)
    [vel_x_aux_dynamic]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
        enable = false  # Disabled initially
    []
    [vel_y_aux_dynamic]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
        enable = false
    []
    # Acceleration computation (dynamic mode only)
    [accel_x_aux]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = 'TIMESTEP_END'
        enable = false
    []
    [accel_y_aux]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = 'TIMESTEP_END'
        enable = false
    []
    [vel_z_zero]
        type = FunctionAux
        variable = vel_z
        function = zero_fn
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [vel_mag_aux]
        type = FarmsVelocityMagnitudeAux
        variable = vel_mag
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
        execute_on = 'TIMESTEP_END'
    []
    # Material property outputs
    [get_xi]
        type = MaterialRealAux
        variable = xi_aux
        property = strain_invariant_ratio
        block = '1 3'
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_aux
        property = second_elastic_strain_invariant
        block = '1 3'
    []
    [get_deviatroic_strain_rate]
        type = MaterialRealAux
        variable = deviatroic_strain_rate_aux
        property = deviatroic_strain_rate
        block = '1 3'
    []
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
    []
    [get_pk2_stress_01]
        type = MaterialRankTwoTensorAux
        variable = pk2_stress_01
        property = pk2_stress
        i = 0
        j = 1
        block = '1 3'
    []
    [get_green_lagrange_elastic_strain_01]
        type = MaterialRankTwoTensorAux
        variable = green_lagrange_elastic_strain_01
        property = green_lagrange_elastic_strain
        i = 0
        j = 1
        block = '1 3'
    []
    [get_plastic_strain_01]
        type = MaterialRankTwoTensorAux
        variable = plastic_strain_01
        property = plastic_strain
        i = 0
        j = 1
        block = '1 3'
    []
    [get_total_lagrange_strain_01]
        type = MaterialRankTwoTensorAux
        variable = total_lagrange_strain_01
        property = total_lagrange_strain
        i = 0
        j = 1
        block = '1 3'
    []
    ###
    [get_eqstrain_nonlocal_initial]
        type = SolutionAux
        variable = eqstrain_nonlocal_initial_aux
        solution = init_sol_components
        from_variable = xi_output
        execute_on = 'INITIAL'
    []
[]

[Kernels]
    # Always active: stress divergence
    [disp_x_kernel]
        type = TotalLagrangianStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
    [disp_y_kernel]
        type = TotalLagrangianStressDivergence
        variable = disp_y
        component = 1
        large_kinematics = true
    []

    # DYNAMIC MODE: Inertia terms (disabled initially, enabled when strain rate > 1e-5)
    # APPLIED TO ALL REGIONS (blocks 1, 2, 3)
    [inertia_x]
        type = InertialForce
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
        enable = false
        block = '1 2 3'
    []
    [inertia_y]
        type = InertialForce
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
        eta = 0
        enable = false
        block = '1 2 3'
    []

    # QUASI-DYNAMIC MODE: Radiation damping
    # Inner region (blocks 1, 3): enabled initially, disabled when strain rate > 1e-5
    [rad_damp_x_inner]
        type = FarmsRadiationDamping
        variable = disp_x
        eta_constant = 1.85e7  # 2 * mu / c_s
        enable = true
        block = '1 2 3'
    []
    [rad_damp_y_inner]
        type = FarmsRadiationDamping
        variable = disp_y
        eta_constant = 1.85e7
        enable = true
        block = '1 2 3'
    []

    # Always active: Rayleigh damping for numerical stability
    [rayleigh_damp_x]
        type = LagrangianStiffPropDampingImplicit
        variable = disp_x
        component = 0
        zeta = 0.04
    []
    [rayleigh_damp_y]
        type = LagrangianStiffPropDampingImplicit
        variable = disp_y
        component = 1
        zeta = 0.04
    []
[]

[Functions]
    [func_top_bc]
        type = ParsedFunction
        expression = 'if (t>dt, 1e-8 * t, 0)'
        symbol_names = 'dt'
        symbol_values = '1e-3'
    []
    #**Parameters:**
    #- Shear modulus: μ = 32.04 GPa
    #    - Plate velocity: V = 1e-8 m/s (typical tectonic rate)
    #- Domain width: L = 60 km = 60,000 m

    #**Result:**
    #```
    #τ̇ = 32.04e9 × 1e-8 / 60000 = 5.34 Pa/s
    #```
    [tectonic_shear_stress]
        type = ParsedFunction
        expression = 'if (t>dt, 13e6 + 5.34 * t, 13e6)'  # Linear stress increase with time
        symbol_names = 'dt'
        symbol_values = '2e-2'
    []
    [zero_fn]
        type = ParsedFunction
        expression = 0
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
    []
    [damage_mat]
        type = DiffusedDamageBreakageMaterialMainApp
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        structural_stress_coefficient = structural_stress_coefficient_aux
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
    []
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Diffused
        large_kinematics = true
        output_properties = 'equiv_plastic_strain equiv_plastic_strain_rate'
        outputs = exodus
        block = '1 3'
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage shear_stress_perturbation damage_perturbation mean_stress_perturbation'
        prop_values = '0.0 0.0 0.0 0.0'
    []
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        block = '2'
    []
    [comp_strain_invariant_ratio]
        type = ComputeXi
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
        block = '2'
    []
    [eqstrain_nonlocal_initial_xi]
      type = CoupledVariableValueMaterial #this material object is in thermalhydraulicApp
      coupled_variable = eqstrain_nonlocal_initial_aux
      prop_name = eqstrain_nonlocal_initial
      output_properties = 'eqstrain_nonlocal_initial'
      outputs = exodus
    []
    [nonlocal_eqstrain]
        type = ElkNonlocalEqstrain
        average_UO = eqstrain_averaging
        output_properties = 'eqstrain_nonlocal'
        outputs = exodus
    []
[]

[UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 200
        prop_name = strain_invariant_ratio
        radius = 400
        weights = BAZANT
        execute_on = TIMESTEP_END
    []
[]

[Preconditioning]
    [smp]
        type = SMP
        full = true
    []
[]

[Postprocessors]
    [_dt]
        type = TimestepSize
    []
    # Velocity monitoring
    [maxvelx]
        type = NodalExtremeValue
        variable = vel_x
    []
    [maxvely]
        type = NodalExtremeValue
        variable = vel_y
    []
    [maxvelz]
        type = NodalExtremeValue
        variable = vel_z
    []
    [maxvelmag]
        type = NodalExtremeValue
        variable = vel_mag
    []
    # KEY: Deviatoric strain rate monitoring for switching
    [max_dev_strain_rate]
        type = ElementExtremeValue
        variable = deviatroic_strain_rate_aux
        value_type = max
    []
    # Adaptive time step bound based on current mode
    [adaptive_dt_bound]
        type = FarmsAdaptiveTimeStepBound
        criterion_postprocessor = max_dev_strain_rate
        threshold = 1e-5
        comparison_type = greater_than
        # When strain rate > 1e-5 (dynamic mode): use small dt_max
        dt_bound_above_threshold = 1e-2
        # When strain rate <= 1e-5 (quasi-dynamic mode): use large dt_max
        dt_bound_below_threshold = 50.0
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    start_time = -1e-12
    end_time = 1e5
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 30
    nl_abs_tol = 1e-10
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
    automatic_scaling = true
    verbose = true
    fixed_point_max_its = 10
    accept_on_max_fixed_point_iteration = false
    fixed_point_rel_tol = 1e-6
    fixed_point_abs_tol = 1e-8

    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-2
        cutback_factor_at_failure = 0.5
        optimal_iterations = 20
        growth_factor = 1.05
        # Use adaptive time step bound from postprocessor
        max_time_step_bound_postprocessor = adaptive_dt_bound
        # Fallback constant value (not used when postprocessor is active)
        max_time_step_bound = 100.0
    []

    [TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    []
[]

[Controls]
    # Disable all dynamic features for the first time step
    [initial_period]
        type = TimePeriod
        start_time = -1e-12
        end_time = 1e-2
        disable_objects = 'AuxKernels/vel_x_aux_dynamic AuxKernels/vel_y_aux_dynamic
                           AuxKernels/accel_x_aux AuxKernels/accel_y_aux
                           Kernels/inertia_x Kernels/inertia_y'
    []

    # BIDIRECTIONAL SWITCHING: Based on deviatoric strain rate
    # When max(strain_rate) > 1e-5: Switch ALL REGIONS to DYNAMIC mode
    # When max(strain_rate) <= 1e-5: Switch ALL REGIONS to QUASI-DYNAMIC mode
    [strain_rate_switch]
        type = FarmsConditionalPostprocessorEnableControl
        postprocessor = max_dev_strain_rate
        comparison_type = greater_than
        threshold = 1e-5
        reverse_on_false = true

        # Enable when strain rate > 1e-5 (DYNAMIC mode in ALL regions)
        enable_objects = 'Kernels/inertia_x Kernels/inertia_y
                          AuxKernels/accel_x_aux AuxKernels/accel_y_aux
                          AuxKernels/vel_x_aux_dynamic AuxKernels/vel_y_aux_dynamic'

        # Disable when strain rate > 1e-5 (turn off QUASI-DYNAMIC in ALL regions)
        disable_objects = 'Kernels/rad_damp_x_inner Kernels/rad_damp_y_inner
                           AuxKernels/vel_x_aux_quasi AuxKernels/vel_y_aux_quasi'
    []
[]

[Outputs]
    [exodus]
        type = Exodus
        time_step_interval = 100
        show = 'disp_x disp_y vel_x vel_y vel_mag alpha_damagedvar_aux B_damagedvar_aux
                xi_aux nonlocal_xi pk2_stress_01 green_lagrange_elastic_strain_01
                plastic_strain_01 total_lagrange_strain_01 deviatroic_strain_rate_aux
                equiv_plastic_strain equiv_plastic_strain_rate'
    []
    #[csv]
    #    type = CSV
    #    time_step_interval = 10
    #    execute_on = 'timestep_end'
    #[]
    [out]
        type = Checkpoint
        time_step_interval = 2000
        num_files = 2
    []
[]

[BCs]
    # Loading boundary conditions
    [preset_displacements]
        type = PresetDisplacement
        boundary = top
        variable = disp_x
        beta = 0.25
        velocity = vel_x
        acceleration = accel_x
        function = func_top_bc
    []
    [initial_shear_stress_top]
        type = NeumannBC
        variable = disp_x
        value = 13e6
        boundary = top
    []
    [static_pressure_top]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -50e6
        displacements = 'disp_x disp_y'
    []
    [static_pressure_bottom]
        type = NeumannBC
        variable = disp_y
        boundary = bottom
        value = 50e6
        displacements = 'disp_x disp_y'
    []
    [static_pressure_left]
        type = NeumannBC
        variable = disp_x
        boundary = left
        value = 50e6
        displacements = 'disp_x disp_y'
    []
    [static_pressure_right]
        type = NeumannBC
        variable = disp_x
        boundary = right
        value = -50e6
        displacements = 'disp_x disp_y'
    []

    # Fixed corner points
    [fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    #add dampers
    [damp_top_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 0
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_top_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 1
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
   [damp_bottom_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 0
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_bottom_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 1
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_left_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 0
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_left_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 1
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_right_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 0
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
    [damp_right_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y'
        velocities = 'vel_x vel_y'
        accelerations = 'accel_x accel_y'
        component = 1
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []
[]

[MultiApps]
    [sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'sub.i'
        execute_on = 'TIMESTEP_END'
        sub_cycling = true
        clone_parent_mesh = true
    []
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_damagedvar_sub B_damagedvar_sub structural_stress_coefficient_sub'
        variable = 'alpha_damagedvar_aux B_damagedvar_aux structural_stress_coefficient_aux'
        execute_on = 'TIMESTEP_END'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'I2_aux nonlocal_xi deviatroic_strain_rate_aux'
        variable = 'I2_sub_aux xi_sub_aux deviatroic_strain_rate_sub_aux'
        execute_on = 'TIMESTEP_END'
    []
[]

[UserObjects]
    [init_sol_components]
        type = SolutionUserObject
        mesh = '../../static_solve/static_solve_large_out.e'
        system_variables = 'disp_x disp_y xi_output I2_output alpha_damagedvar_output B_damagedvar_output'
        timestep = LATEST
        force_preaux = true
        execute_on = 'INITIAL'
    []
[]

[ICs]
    [disp_x_ic]
        type = SolutionIC
        variable = disp_x
        solution_uo = init_sol_components
        from_variable = disp_x
    []
    [disp_y_ic]
        type = SolutionIC
        variable = disp_y
        solution_uo = init_sol_components
        from_variable = disp_y
    []
    [strain_invariant_ratio_ic]
        type = SolutionIC
        variable = nonlocal_xi
        solution_uo = init_sol_components
        from_variable = xi_output
    []
    [I2_aux_ic]
        type = SolutionIC
        variable = I2_aux
        solution_uo = init_sol_components
        from_variable = I2_output
    []
    [alpha_damagedvar_sub_ic]
        type = SolutionIC
        variable = alpha_damagedvar_aux
        solution_uo = init_sol_components
        from_variable = alpha_damagedvar_output
    []
    [B_damagedvar_sub_ic]
        type = SolutionIC
        variable = B_damagedvar_aux
        solution_uo = init_sol_components
        from_variable = B_damagedvar_output
    []
[]

[VectorPostprocessors]
  [top_surface_vel]
    type = LineValueSampler
    start_point = '-15000 100 0'
    end_point = '15000 100 0'
    sort_by = x
    num_points = 300
    variable = 'vel_x vel_y disp_x disp_y'
  []
  [bottom_surface_vel]
    type = LineValueSampler
    start_point = '-15000 -100 0'
    end_point = '15000 -100 0'
    sort_by = x
    num_points = 300
    variable = 'vel_x vel_y disp_x disp_y'
  []
[]
