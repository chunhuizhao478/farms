fluid_density = 1000   
solid_density = 2700
gravity_pos = 9.81

#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/mesh_large.msh'
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
        # coord = ' -120000 -120000 -120000;
        #            120000 -120000 -120000;
        #            120000 120000  -120000;
        #           -120000 120000  -120000'
        coord = ' -120000 -120000 -120000;
        120000 120000  -120000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9
        
    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.7
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.8
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
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
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
    #
    [alpha_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    [B_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
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
    #
    [gradx_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    [grady_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    #spatial damage parameters
    [cg_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
    [nonlocal_xi]
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
    #
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
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
    #
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
[]

#parameters for the initial stress field
################################################
bxx = 0.926793
byy = 1.073206
bxy = -0.169029
linear_variation_cutoff_distance = 15600
################################################
[Functions]
    [func_pos_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    #
    [func_pos_zz_stress]
        type = InitialDepthDependentStress
        i = 3
        j = 3
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_zz_stress]
        type = InitialDepthDependentStress
        i = 3
        j = 3
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
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
        # outputs = exodus
    []
    # # damage
    [damage_mat]
        type = DiffusedDamageBreakageMaterialMainApp
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        structural_stress_coefficient = structural_stress_coefficient_aux
        #build L matrix using velocity
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Diffused
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain strain_invariant_ratio'
        outputs = exodus
        block = '1 3'
    []
    [dummy_initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = '0.0'
    []
    #elastic material
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
        block = '2'
    []
    #strain invariant ratio
    [comp_strain_invariant_ratio]
        type = ComputeXi 
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
        block = '2'
    []
    #nonlocal eqstrain
    [nonlocal_eqstrain]
        type = ElkNonlocalEqstrain
        average_UO = eqstrain_averaging
        output_properties = 'eqstrain_nonlocal'
        outputs = exodus
    []
    #shear stress perturbation
    [damage_perturbation]
        type = PerturbationRadial
        nucl_center = '0 0 0'
        peak_value = 0
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'shear_stress'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
        outputs = exodus
    []
[] 

[UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 300
        prop_name = strain_invariant_ratio
        radius = 200
        weights = BAZANT
        execute_on = LINEAR
    []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */vel_z */accel_x */accel_y */accel_z */inertia_x */inertia_y */inertia_z */damp_left_x */damp_left_y */damp_left_z */damp_right_x */damp_right_y */damp_right_z */damp_bottom_x */damp_bottom_y */damp_bottom_z */damp_front_x */damp_front_y */damp_front_z */damp_back_x */damp_back_y */damp_back_z'
    start_time = -1e-12
    end_time = 1e-2 # dt used in the simulation
  []
[../]
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 1e10
    num_steps = 50
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 10
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
        dt = 1e-2
        cutback_factor_at_failure = 0.5
        optimal_iterations = 8
        growth_factor = 1.1
        max_time_step_bound = 1e7
        #constrain velocity during dynamic simulation
        constrain_by_velocity = true
        vel_threshold = 1e-2
        constant_dt_on_overspeed = 1e-2
        maxvelx = 'maxvelx'
        maxvely = 'maxvely'
    []
    # [TimeStepper]
    #     type = IterationAdaptiveDT
    #     cutback_factor_at_failure = 0.5
    #     growth_factor = 2.0
    #     optimal_iterations = 100
    #     dt = 1e-2
    # []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

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
        time_step_interval = 10
        show = 'vel_x vel_y alpha_damagedvar_aux B_damagedvar_aux xi_aux deviatroic_strain_rate_aux nonlocal_xi pk2_stress_01 green_lagrange_elastic_strain_01 plastic_strain_01 total_lagrange_strain_01'
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
    [../]
[]

#apply boundary conditions
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

[BCs]
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
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
        shear_wave_speed = 3464
        p_wave_speed = 6000
        density = 2700
    []       
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'dynamic_solve_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
        sub_cycling = true
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_damagedvar_sub B_damagedvar_sub structural_stress_coefficient_sub'
        variable = 'alpha_damagedvar_aux B_damagedvar_aux structural_stress_coefficient_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'I2_aux nonlocal_xi deviatroic_strain_rate_aux'
        variable = 'I2_sub_aux xi_sub_aux deviatroic_strain_rate_sub_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_test_cdbm_out.e'
      system_variables = 'disp_x disp_y disp_z initial_xi_aux initial_I2_aux initial_damage_aux initial_breakage_aux'
      timestep = LATEST
      force_preaux = true
      execute_on = 'INITIAL'
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
    [strain_invariant_ratio_ic]
      type = SolutionIC
      variable = nonlocal_xi
      solution_uo = init_sol_components
      from_variable = initial_xi_aux
    []
    [I2_aux_ic]
      type = SolutionIC
      variable = I2_aux
      solution_uo = init_sol_components
      from_variable = initial_I2_aux
    []
    [alpha_damagedvar_sub_ic]
        type = SolutionIC
        variable = alpha_damagedvar_aux
        solution_uo = init_sol_components
        from_variable = initial_damage_aux
    []  
    [B_damagedvar_sub_ic]
        type = SolutionIC
        variable = B_damagedvar_aux
        solution_uo = init_sol_components
        from_variable = initial_breakage_aux
    []    
[]