[Mesh]
    [./msh]
        type = FileMeshGenerator
        # file = '../meshfile/tpv2052dm.msh'
        file = '../../meshfile/tpv2052dm_2ndorder_mirrormesh_small_local.msh'
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
        coord = '0 -7500 0'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y'
[]

[GlobalParams]

    displacements = 'disp_x disp_y'

    use_displaced_mesh = false
    
    ##----continuum damage breakage model----##
    #initial lambda value (SECOND lame constant) [Pa]
    lambda_o = 10e9
        
    #initial shear modulus value (SECOND lame constant) [Pa]
    shear_modulus_o = 10e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.9
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 300

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 500

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 10

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 3

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.01 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.7
    
[]


[Variables]
    [disp_x]
        order = SECOND
        family = LAGRANGE     
    []
    [disp_y]
        order = SECOND
        family = LAGRANGE    
    []
[]

[AuxVariables]
    [vel_x]
        order = SECOND
        family = LAGRANGE
    []
    [accel_x]
        order = SECOND
        family = LAGRANGE
    []
    [vel_y]
        order = SECOND
        family = LAGRANGE
    []
    [accel_y]
        order = SECOND
        family = LAGRANGE
    []
    #
    [alpha_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [B_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [strain_invariant_ratio_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [correlated_xio]
        order = SECOND
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
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    #
    [alpha_damagedvar_aux]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
        execute_on = 'timestep_end'
        block = '1 2'
    []
    [B_damagedvar_aux]
        type = MaterialRealAux
        variable = B_damagedvar_aux
        property = B_damagedvar
        execute_on = 'timestep_end'
        block = '1 2'
    []  
    [strain_invariant_ratio_aux]
        type = MaterialRealAux
        variable = strain_invariant_ratio_aux
        property = strain_invariant_ratio
        execute_on = 'timestep_end'
        block = '1 2'
    []
    [get_initial_damage]
        type = SolutionAux
        variable = initial_damage_aux
        solution = init_sol_components
        from_variable = initial_damage_aux
    []
    #
    [get_correlated_xio]
        type = FunctionAux
        variable = correlated_xio
        function = node
        execute_on = 'INITIAL TIMESTEP_END'
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
    # # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar xi_0'
        outputs = exodus
        block = '1 2'
        # Option 1: Use constant value
        # use_xi0_aux = false 
        # Option 2: Use aux variable
        use_xi0_aux = true
        xi0_aux = correlated_xio
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2
        large_kinematics = true
        output_properties = 'pk1_stress pk2_stress green_lagrange_strain plastic_strain deviatroic_stress strain_invariant_ratio'
        outputs = exodus
        block = '1 2'
    []
    # elastic
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
        block = 3
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
        block = 3
    []
    [define_initial_damage_matprop]
        type = ParsedMaterial
        property_name = initial_damage
        coupled_variables = 'initial_damage_aux'
        expression = 'initial_damage_aux'
        outputs = exodus
    []
[]  

[Functions]
    [func_top_bc]
        type = ParsedFunction
        expression = 'if (t>dt, 1e-8 * t, 0)'
        symbol_names = 'dt'
        symbol_values = '1e-2'
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
    num_steps = 1
    l_max_its = 100
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
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    # dt = 1e-2
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 0.01
        cutback_factor_at_failure = 0.5
        optimal_iterations = 8
        growth_factor = 1.5
        max_time_step_bound = 1e10
    []
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
    disable_objects = '*/vel_x */vel_y */accel_x */accel_y */inertia_x */inertia_y */bc_load_top_x */damp_left_x */damp_left_y */damp_right_x */damp_right_y'
    start_time = -1e-12
    end_time = 1e-2 # dt used in the simulation
  []
[../]

[Postprocessors]
    [./_dt]
        type = TimestepSize
    [../]
    [./maxvelx]
        type = NodalExtremeValue
        variable = vel_x
    [../]
    [./maxvely]
        type = NodalExtremeValue
        variable = vel_y
    [../]
[../]

[Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 1
      show = 'vel_x vel_y initial_damage alpha_damagedvar_aux B_damagedvar_aux strain_invariant_ratio_aux pk2_stress_00 pk2_stress_11 pk2_stress_01 correlated_xio xi_0'
    [../]
    [./csv]
      type = CSV
      time_step_interval = 1
      show = '_dt maxvelx maxvely'
    [../]
[]

[BCs]
    [bc_load_top_x]
        type = PresetDisplacement
        boundary = top
        variable = disp_x
        beta = 0.25
        velocity = vel_x
        acceleration = accel_x
        function = func_top_bc
    []
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    [] 
    [./Pressure]
        [static_pressure_top]
            boundary = top
            factor = 120e6
            displacements = 'disp_x disp_y'
        []    
        [static_pressure_left]
            boundary = left
            factor = 135e6
            displacements = 'disp_x disp_y'
        []  
        [static_pressure_right]
            boundary = right
            factor = 135e6
            displacements = 'disp_x disp_y'
        []     
    []        
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr2_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    #add dampers
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
        shear_wave_speed = 1924.5
        p_wave_speed = 3333.3
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
        shear_wave_speed = 1924.5
        p_wave_speed = 3333.3
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
        shear_wave_speed = 1924.5
        p_wave_speed = 3333.3
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
        shear_wave_speed = 1924.5
        p_wave_speed = 3333.3
        density = 2700
    []
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = './static_solve_out.e'
      system_variables = 'disp_x disp_y initial_damage_aux'
      timestep = LATEST
      force_preaux = true
    [../]
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
[]

[UserObjects]
    [reader_node]
        type = PropertyReadFile
        prop_file_name = '../mapped_weibull_field.csv'
        read_type = 'node'
        nprop = 3 # number of columns in CSV
    []
[]

[Functions]
    [node]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'reader_node'
        read_type = 'node'
        # 0-based indexing
        column_number = '2'
    []
[]