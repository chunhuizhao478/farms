#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/mesh_local.msh'
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
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9
        
    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9
    
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
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8

    #add strain rate dependent Cd option
    m_exponent = 0.8
    strain_rate_hat = 1e-8
    cd_hat = 1e3
    
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
    [vel_z]
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
    #
    [initial_damage_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [deviatroic_strain_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [Cd_rate_dependent_aux]
        order = FIRST
        family = MONOMIAL
    []
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
    #
    [nonlocal_xi]
        order = FIRST
        family = MONOMIAL
    []
[]

[AuxKernels]
    #options to use Newmark Beta method
    #--------------------------------------------------------------#
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
    #--------------------------------------------------------------#
    #output material properties
    [get_deviatroic_strain_rate]
        type = MaterialRealAux
        variable = deviatroic_strain_rate_aux
        property = deviatroic_strain_rate
        block = '1 3'
    []
    [get_cd_rate_dependent]
        type = MaterialRealAux
        variable = Cd_rate_dependent_aux
        property = Cd_rate_dependent
        block = '1 3'
    []
    [get_alpha_damagedvar]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
        block = '1 3'
    []
    [get_B_damagedvar]
        type = MaterialRealAux
        variable = B_damagedvar_aux
        property = B_damagedvar
        block = '1 3'
    []
    [get_strain_invariant_ratio]
        type = MaterialRealAux
        variable = strain_invariant_ratio_aux
        property = strain_invariant_ratio
        block = '1 3'
    []
    #
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
    []
    #--------------------------------------------------------------#
    #aux parameters for damage breakage model
    [get_cd_block13]
        type = ConstantAux
        variable = Cd_constant_aux
        value = 0
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
    #--------------------------------------------------------------#
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
        outputs = exodus
    []
    # # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar'
        outputs = exodus
        #options to use auxiliary variables
        #-------------------------------------------------------#
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
        #-------------------------------------------------------#
        #options to use strain rate dependent Cd
        use_cd_strain_dependent = true
        #m_exponent = 0.8
        #strain_rate_hat = 1e-8
        #cd_hat = 1e3
        #options to use velocity to build L matrix
        use_vels_build_L = true
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
        #--------------------------------------------------------#
        #options to use nonlocal xi
        use_nonlocal_xi = true
        nonlocal_xi = nonlocal_xi
        #-------------------------------------------------------#
        # use initial damage time dependent
        build_param_use_initial_damage_time_dependent_mat = true
        build_param_peak_value = 0.7
        build_param_sigma = 5e2
        build_param_len_of_fault = 8000
        #-------------------------------------------------------#
        # use damage perturbation time dependent
        perturbation_build_param_use_damage_perturb = true
        perturbation_build_param_nucl_center = '0 0'
        perturbation_build_param_length = 8000
        perturbation_build_param_thickness = 200
        perturbation_build_param_peak_value = 0.3
        perturbation_build_param_sigma = 1318.02
        perturbation_build_param_duration = 1.0
        #-------------------------------------------------------#
        block = '1 3'
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Debug
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain strain_invariant_ratio'
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
[] 

[UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 200
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
  
[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = 30
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 5
    nl_abs_tol = 1e-10
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

[Controls] # turns off inertial terms for the FIRST time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */accel_x */accel_y */inertia_x */inertia_y'
    start_time = -1e-12
    end_time = 1e-2 # dt used in the simulation
  []
[../]

[Outputs]
    #save the solution to a exodus file every 0.1 seconds
    [./exodus]
      type = Exodus
      time_step_interval = 5
      show = 'vel_x vel_y alpha_damagedvar_aux B_damagedvar_aux strain_invariant_ratio_aux pk2_stress_01 green_lagrange_elastic_strain_01 plastic_strain_01 deviatroic_strain_rate_aux Cd_rate_dependent_aux nonlocal_xi' 
    [../]
    # #save the solution to a csv file every 0.001 seconds
    # [./csv]
    #   type = CSV
    #   time_step_interval = 1
    # [../]
[]

[BCs]
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    #add initial shear stress
    [./initial_shear_stress]
        type = NeumannBC
        variable = disp_x
        value = 12e6
        boundary = top
    [] 
    # 
    [static_pressure_top]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -50e6
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
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_test_out.e'
      system_variables = 'disp_x disp_y'
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