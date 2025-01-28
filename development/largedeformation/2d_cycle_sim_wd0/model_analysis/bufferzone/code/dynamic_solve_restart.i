#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/mesh.msh'
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
        coord = '0 -30000 0'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y'
[]

[GlobalParams]

    displacements = 'disp_x disp_y'
    
    ##----continuum damage breakage model----##
    #initial lambda value (SECOND lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (SECOND lame constant) [Pa]
    shear_modulus_o = 30e9
    
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

    #if option 2, use Cd_constant #specify by auxiliary variable
    # Cd_constant = 10

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd #specify by auxiliary variable
    # CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0 #specify by auxiliary variable
    # CBH_constant = 10

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf" #specify by auxiliary variable
    # C_1 = 3

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    # C_2 = 0.05

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

    #add strain rate dependent Cd option
    # m_exponent = 0.85
    # strain_rate_hat = 1e-4
    # cd_hat = 300 #incresed from 1 to 100
    
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
    # [alpha_damagedvar_aux_SECONDmono]
    #     order = SECOND
    #     family = MONOMIAL
    # []
    # [B_damagedvar_aux_SECONDmono]
    #     order = SECOND
    #     family = MONOMIAL
    # []
    # [strain_invariant_ratio_aux_SECONDmono]
    #     order = SECOND
    #     family = MONOMIAL
    # []
    # [initial_damage_aux_SECONDmono]
    #     order = SECOND
    #     family = MONOMIAL 
    # []
    #
    # [timeintegratorflag]
    #     order = SECOND
    #     family = MONOMIAL
    # []
    #
    [Cd_constant_aux]
        order = CONSTANT
        family = MONOMIAL
    []  
    [Cb_multiplier_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    [Cbh_constant_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    [C1_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    [C2_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [strain_invariant_ratio_const_aux]
        order = CONSTANT
        family = MONOMIAL
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
    #get damage,breakage,strain invariant ratio in constant monomial
    [alpha_damagedvar_aux]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
        execute_on = 'timestep_end'
    []
    [B_damagedvar_aux]
        type = MaterialRealAux
        variable = B_damagedvar_aux
        property = B_damagedvar
        execute_on = 'timestep_end'
    []  
    [strain_invariant_ratio_aux]
        type = MaterialRealAux
        variable = strain_invariant_ratio_aux
        property = strain_invariant_ratio
        execute_on = 'timestep_end'
    []
    # #get damage,breakage,strain invariant ratio in SECOND monomial
    # [alpha_damagedvar_aux_SECONDmono]
    #     type = MaterialRealAux
    #     variable = alpha_damagedvar_aux_SECONDmono
    #     property = alpha_damagedvar
    #     execute_on = 'timestep_end'
    # []
    # [B_damagedvar_aux_SECONDmono]
    #     type = MaterialRealAux
    #     variable = B_damagedvar_aux_SECONDmono
    #     property = B_damagedvar
    #     execute_on = 'timestep_end'
    # []  
    # [strain_invariant_ratio_aux_SECONDmono]
    #     type = MaterialRealAux
    #     variable = strain_invariant_ratio_aux_SECONDmono
    #     property = strain_invariant_ratio
    #     execute_on = 'timestep_end'
    # []    
    # #
    # [get_flag]
    #     type = MaterialRealAux
    #     variable = timeintegratorflag
    #     property = flag
    #     execute_on = 'timestep_end'
    # []
    #block 1: inner block where damage is activated
    #block 2: outer block where damage is not activated
    #block 3: buffer block where damage/breakage accumulation is decreased, using xi only
    #cd constant
    [get_cd_block13]
        type = ConstantAux
        variable = Cd_constant_aux
        value = 300
        block = '1 3'
        execute_on = 'INITIAL'
    []
    [get_cd_block23]
        type = ConstantAux
        variable = Cd_constant_aux
        value = 0
        block = '2'
        execute_on = 'INITIAL'
    []    
    #cb multiplier
    [get_cb_multiplier_block13]
        type = ConstantAux
        variable = Cb_multiplier_aux
        value = 500
        block = '1 3'
        execute_on = 'INITIAL'
    []
    [get_cb_multiplier_block2]
        type = ConstantAux
        variable = Cb_multiplier_aux
        value = 0
        block = '2'
        execute_on = 'INITIAL'
    []
    #cbh constant 
    [get_cbh_constant_block1]
        type = ConstantAux
        variable = Cbh_constant_aux
        value = 1e4
        block = '1 3'
        execute_on = 'INITIAL'
    []
    [get_cbh_constant_block2]
        type = ConstantAux
        variable = Cbh_constant_aux
        value = 0
        block = 2
        execute_on = 'INITIAL'
    []
    #C1
    [get_c1_block13]
        type = ConstantAux
        variable = C1_aux
        value = 300
        block = '1 3'
        execute_on = 'INITIAL'
    []
    [get_c1_block2]
        type = ConstantAux
        variable = C1_aux
        value = 0
        block = 2
        execute_on = 'INITIAL'
    []
    #C2
    [get_c2_block13]
        type = ConstantAux
        variable = C2_aux
        value = 0.05
        block = '1 3'
        execute_on = 'INITIAL'
    []
    [get_c2_block2]
        type = ConstantAux
        variable = C2_aux
        value = 0
        block = '2'
        execute_on = 'INITIAL'
    []
    [get_initial_damage]
        type = SolutionAux
        variable = initial_damage_aux
        solution = init_sol_components
        from_variable = initial_damage_aux
    []
    #strain invariant ratio xi
    #the constant strain invariant ratio xi is only activated in the buffer block
    #in the material object, the block where this constant xi is activated must be speficied
    [get_xi_block3]
        type = ConstantAux
        variable = strain_invariant_ratio_const_aux
        value = -1.7
        block = '3'
        execute_on = 'INITIAL'
    []
    #don't apply this in the material object, as block 1
    #xi needs to be actually calculated
    [get_xi_block12]
        type = ConstantAux
        variable = strain_invariant_ratio_const_aux
        value = 0
        block = '1 2'
        execute_on = 'INITIAL'
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
        zeta = 0.5
    []
    [damping_y]
        type = StiffPropDampingImplicit
        variable = disp_y
        component = 1
        zeta = 0.5
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
        output_properties = 'alpha_damagedvar B_damagedvar'
        outputs = exodus
        #options to use auxiliary variables
        use_cd_aux = true
        Cd_constant_aux = Cd_constant_aux
        use_cb_multiplier_aux = true
        Cb_multiplier_aux = Cb_multiplier_aux
        use_cbh_aux = true
        CBH_aux = Cbh_constant_aux
        use_c1_aux = true
        C1_aux = C1_aux
        use_c2_aux = true
        C2_aux = C2_aux
        # use_cd_strain_dependent = true
        # use_total_strain_rate = true
        # block_id_applied = 1
        use_const_xi_aux = true
        const_xi_aux = strain_invariant_ratio_const_aux
        const_xi_block_id = 3
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2
        large_kinematics = true
        output_properties = 'pk1_stress pk2_stress green_lagrange_elastic_strain plastic_strain deviatroic_stress strain_invariant_ratio total_lagrange_strain Cd_constant_aux'
        outputs = exodus
        # block = '2'
    []
    [define_initial_damage_matprop]
        type = ParsedMaterial
        property_name = initial_damage
        coupled_variables = 'initial_damage_aux'
        expression = 'initial_damage_aux'
        outputs = exodus
    []
    #
    # [forcedampingflag]
    #     type = ForceDampingFlag
    #     vel_maximum_threshold = 1e-2
    #     vel_minimum_threshold = 1e-6
    #     max_vel_x = maxvelx
    #     max_vel_y = maxvely
    # []
[]  

[Functions]
    [func_top_bc]
        type = ParsedFunction
        expression = 'if (t>dt, 1e-8 * t, 0)'
        symbol_names = 'dt'
        symbol_values = '1e-2'
    []
    [func_bottom_bc]
        type = ParsedFunction
        expression = 'if (t>dt, -1e-8 * t, 0)'
        symbol_names = 'dt'
        symbol_values = '1e-2'
    []
    [func_cd]
        type = SpatialDamageBreakageParameters
        xmin = -15000
        xmax = 15000
        ymin = -3000
        ymax = 3000
        max_val = 300
        min_val = 0.1
        scale = 500
    []
    [func_cb_multiplier]
        type = SpatialDamageBreakageParameters
        xmin = -15000
        xmax = 15000
        ymin = -3000
        ymax = 3000
        max_val = 500
        min_val = 0.1
        scale = 500        
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
    start_time = 8376319166.877851
    end_time = 1e100
    # num_steps = 10
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
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 0.01
        cutback_factor_at_failure = 0.5
        optimal_iterations = 8
        growth_factor = 1.5
        max_time_step_bound = 1e10
    []
    # [./TimeIntegrator]
    #     type = FarmsNewmarkBeta
    #     beta = 0.25
    #     gamma = 0.5
    #     factor = 0.9
    #     threshold = 1e-3
    # [../]
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

# [Controls] # turns off inertial terms for the SECOND time step
#   [./period0]
#     type = TimePeriod
#     disable_objects = '*/vel_x */vel_y */accel_x */accel_y */inertia_x */inertia_y */bc_load_top_x */damp_left_x */damp_left_y */damp_right_x */damp_right_y'
#     start_time = -1e-12
#     end_time = 1e-2 # dt used in the simulation
#   []
# [../]

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
    #must be named "flag", this will be called in FarmsNewmarkBeta
    # [flag]
    #     type = ElementExtremeValue
    #     variable = timeintegratorflag
    #     force_postaux = true
    # []
[../]

[Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 100
    #   show = 'disp_x disp_y vel_x vel_y initial_damage alpha_damagedvar_aux B_damagedvar_aux strain_invariant_ratio_aux pk2_stress_00 pk2_stress_11 pk2_stress_01 pk2_stress_22 plastic_strain_00 plastic_strain_01 plastic_strain_11 plastic_strain_22 green_lagrange_elastic_strain_00 green_lagrange_elastic_strain_01 green_lagrange_elastic_strain_11 green_lagrange_elastic_strain_22 deviatroic_stress_00 deviatroic_stress_01 deviatroic_stress_11 deviatroic_stress_22 strain_invariant_ratio total_lagrange_strain_00 total_lagrange_strain_01 total_lagrange_strain_11 total_lagrange_strain_22 Cd_rate_dependent_aux strain_dir0_positive_aux Cd_constant_aux'
    [../]
    [./csv]
      type = CSV
      time_step_interval = 1
      show = '_dt maxvelx maxvely'
    [../]
    [out]
        type = Checkpoint
        time_step_interval = 100
        num_files = 2
    []
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
    # 
    [./Pressure]
        [static_pressure_top]
            boundary = top
            factor = 50e6
            displacements = 'disp_x disp_y'
        []    
        [static_pressure_left]
            boundary = left
            factor = 50e6
            displacements = 'disp_x disp_y'
        []  
        [static_pressure_right]
            boundary = right
            factor = 50e6
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
        shear_wave_speed = 3333
        p_wave_speed = 5773
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
        shear_wave_speed = 3333
        p_wave_speed = 5773
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
        shear_wave_speed = 3333
        p_wave_speed = 5773
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
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
[]

# [UserObjects]
#     [./init_sol_components]
#       type = SolutionUserObject
#       mesh = './static_solve_out.e'
#       system_variables = 'initial_damage_aux disp_x disp_y'
#       timestep = LATEST
#       force_preaux = true
#     [../]
# []

# [ICs]
#     [disp_x_ic]
#       type = SolutionIC
#       variable = disp_x
#       solution_uo = init_sol_components
#       from_variable = disp_x
#     []
#     [disp_y_ic]
#       type = SolutionIC
#       variable = disp_y
#       solution_uo = init_sol_components
#       from_variable = disp_y
#     []
# []

[Problem]
    #Note that the suffix is left off in the parameter below.
    restart_file_base = dynamic_solve_out_cp/LATEST  # You may also use a specific number here
[]