#continuum damage-breakage model dynamics

##########################################################################################################################################
#Mesh section
#FileMeshGenerator: read mesh file
#SideSetsFromNormalsGenerator: generate side sets from normals
#ExtraNodesetGenerator: generate extra nodeset - here we use it to define corner points associated with the bottom boundary
##########################################################################################################################################
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_large_extended.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0
                    0 0 -1
                    0 0 1'
        new_boundary = 'left right front back bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = ' -12000 -10000 -20000;
                   12000 -10000 -20000;
                   12000 10000  -20000;
                  -12000 10000  -20000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.8
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    # #if option 2, use Cd_constant
    # Cd_constant = 1e4

    # #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    # #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    # CdCb_multiplier = 1000

    # #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # # CBCBH_multiplier = 0.0
    # CBH_constant = 1e4

    # #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    # C_1 = 300

    # #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    # C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    # energy ratio
    chi = 0.7
    
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
    [disp_z]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
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
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
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
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [alpha_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [B_aux]
        order = FIRST
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
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    #
    #aux parameters for damage breakage model
    [get_cd_block13]
        type = ConstantAux
        variable = Cd_constant_aux
        value = 1e4
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
        value = 1000
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
    #
    [get_initial_damage]
        type = MaterialRealAux
        variable = initial_damage_aux
        property = initial_damage_time_dependent_material
    []
    [get_damage]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
        block = '1 3'
    []
    [get_strain_invariant_ratio]
        type = MaterialRealAux
        variable = xi_aux
        property = strain_invariant_ratio
        block = '1 3'
    []
    [get_B]
        type = MaterialRealAux
        variable = B_aux
        property = B_damagedvar
        block = '1 3'
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
    [dispkernel_z]
        type = TotalLagrangianStressDivergence
        variable = disp_z
        component = 2
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
    [./inertia_z]
        type = InertialForce
        variable = disp_z
        acceleration = accel_z
        velocity = vel_z
        beta = 0.25
        gamma = 0.5
        eta = 0
    []  
    [gravity]
        type = Gravity
        variable = disp_z
        value = -9.81
    []       
[]

[Materials]
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        outputs = exodus
    []
    [density]
        type = GenericConstantMaterial
        prop_names = 'density nonADdensity'
        prop_values = '2700 2700'
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
        # use initial damage time dependent
        build_param_use_initial_damage_time_dependent_mat = true
        build_param_peak_value = 0.7
        build_param_sigma = 5e2
        build_param_len_of_fault = 14000
        build_param_use_initial_damage_3D = true
        build_param_len_of_fault_dip = 10000
        build_param_center_point = '0 0 -10000'
        # use damage perturbation time dependent
        perturbation_build_param_use_damage_perturb = true
        perturbation_build_param_nucl_center = '0 0 -7500'
        perturbation_build_param_length = 1000
        perturbation_build_param_thickness = 200
        perturbation_build_param_peak_value = 0.3
        perturbation_build_param_sigma = 2.0
        perturbation_build_param_duration = 0.01
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Debug
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain strain_invariant_ratio'
        outputs = exodus
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
    []
[]  

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 100.0
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 5
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    dt = 1e-3
    verbose = true
    # [TimeStepper]
    #     type = FarmsIterationAdaptiveDT
    #     dt = 1e-2
    #     cutback_factor_at_failure = 0.5
    #     optimal_iterations = 8
    #     growth_factor = 1.5
    #     max_time_step_bound = 1e10
    # []
    # [./TimeStepper]
    #     type = SolutionTimeAdaptiveDT
    #     dt = 0.01
    # [../]
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

[Controls] # turns off inertial terms for the FIRST time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */vel_z */accel_x */accel_y */accel_z */inertia_x */inertia_y */inertia_z */damp_top_x */damp_top_y */damp_top_z */damp_bottom_x */damp_bottom_y */damp_bottom_z */damp_left_x */damp_left_y */damp_left_z */damp_right_x */damp_right_y */damp_right_z */damp_front_x */damp_front_y */damp_front_z */damp_back_x */damp_back_y */damp_back_z'
    start_time = -1e-12
    end_time = 1e-3 # dt used in the simulation
  []
[../]

[Outputs] 
    ### save the solution to a exodus file every [time_step_interval] time steps]
    exodus = true
    time_step_interval = 100
    show = 'strain_invariant_ratio vel_x vel_y vel_z alpha_damagedvar_aux B_aux xi_aux'
    # [./csv]
    #     type = CSV
    #     time_step_interval = 1
    #     show = 'maxvelx maxvely maxvelz'
    # [../]
[]

[BCs]
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []
    #Note: use neuamnnBC gives minimum waves than pressureBC  
    [static_pressure_left]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = left
        function = func_pos_xx_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = right
        function = func_neg_xx_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    #
    [static_pressure_front]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = front
        function = func_pos_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = back
        function = func_neg_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []
    #
    [static_pressure_front_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = front
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_left_shear]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = left
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right_shear]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = right
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []   
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = DirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]

[Functions]
    [func_pos_yy_stress]
        type = ParsedFunction      
        expression = 'if(-z<15600, -1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_yy_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'  
    []
    [func_pos_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600, -1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'
    []
    [func_pos_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, -1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, -1 * (-0.8 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
    [func_neg_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, 1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, 1 * (-0.8 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
[]

# #####################################
# #BCs Section
# #Absorbing boundary conditions
# #####################################
[BCs]
    ##non-reflecting bc
    #
    #add dampers
    [damp_top_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_top_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_top_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    #
    [damp_bottom_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_bottom_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_bottom_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    #
    [damp_front_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_front_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_front_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    #
    [damp_back_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_back_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_back_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    #
    [damp_left_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
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
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_left_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    #  
    [damp_right_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
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
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
    [damp_right_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333
        p_wave_speed = 5773
        density = 2700
    []
[] 

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_extended_implicit_out.e'
      system_variables = 'disp_x disp_y disp_z'
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
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]
