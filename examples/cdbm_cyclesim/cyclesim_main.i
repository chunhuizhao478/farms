#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  './meshfile/cyclesim.msh'
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
    [./elasticblock_1]
        type = SubdomainBoundingBoxGenerator    
        input = sidesets
        block_id = 1
        bottom_left = '-20 -50 0'
        top_right = '20 -45 0'
    []
    [./elasticblock_2]
        type = SubdomainBoundingBoxGenerator    
        input = elasticblock_1
        block_id = 1
        bottom_left = '-20 45 0'
        top_right = '20 50 0'
    []      
[]

[GlobalParams]
  
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y'
  
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
      
    #initial shear modulus value (second lame constant) [Pa]
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
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
  
    ##Compute gamma_damaged_r, xi_1
    #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
    #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
    #check struct_param.m 
  
    #coefficient of damage solid modulus
    gamma_damaged_r = 34.785e9
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.825
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.8
    a0 = 7.42e9
    a1 = -21.341e9
    a2 = 19.028e9
    a3 = -4.924e9
  
    #diffusion coefficient #for structural stress coupling
    D = 0
  
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
    #damage-breakage model
    [./B_old]
        order = FIRST
        family = LAGRANGE
    []
    [./xi_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./I2_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./lambda_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./gamma_old]
        order = CONSTANT
        family = MONOMIAL
    []
    #updated alpha, B
    [./alpha_in]
        order = CONSTANT
        family = MONOMIAL
    []
    [./B_in]
        order = CONSTANT
        family = MONOMIAL
    []
    #output alpha, B to subApp
    [./alpha_damagedvar_out]
        order = CONSTANT
        family = MONOMIAL
    []
    [./B_out]
        order = CONSTANT
        family = MONOMIAL   
    []        
    #grad_alpha
    [./alpha_grad_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [./alpha_grad_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [./alpha_grad_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #initial alpha
    [./initial_alpha]
        order = CONSTANT
        family = MONOMIAL        
    []
    #principal strain
    [./principal_strain_rate_out]
        order = CONSTANT
        family = MONOMIAL         
    []
    #record shear stress
    [./shear_stress_applied]
        order = CONSTANT
        family = MONOMIAL        
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        zeta = 0.0029
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        zeta = 0.0029
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
        eta = 31.35
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
        eta = 31.35
    []
[]

[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = timestep_end
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = timestep_end
    []
    #
    [get_xi_old]
        type = ADMaterialRealAux
        property = xi
        variable = xi_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_I2_old]
        type = ADMaterialRealAux
        property = I2
        variable = I2_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_mu_old]
        type = ADMaterialRealAux
        property = shear_modulus
        variable = mu_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_lambda_old]
        type = ADMaterialRealAux
        property = lambda
        variable = lambda_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_gamma_old]
        type = ADMaterialRealAux
        property = gamma_damaged
        variable = gamma_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_alpha_old]
        type = ADMaterialRealAux
        property = alpha_damagedvar
        variable = alpha_damagedvar_out
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    [get_B_old]
        type = ADMaterialRealAux
        property = B
        variable = B_out
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    #add inital distribution of alpha
    [assign_initial_alpha]
        type = FunctionAux
        variable = initial_alpha
        function = func_stress_yy
        execute_on = 'INITIAL TIMESTEP_BEGIN'
        block = 0
    []
    #principal strain
    [get_principal_strain]
        type = ADMaterialRateRealAux
        property = principal_strain
        variable = principal_strain_rate_out
        execute_on = 'INITIAL TIMESTEP_END'
        block = 0
    []
    [record_applied_shear_stress]
        type = FunctionAux
        function = func_stress_xy
        variable = shear_stress_applied
    []
[]

[Materials]
    [damagestress]
        type = ADComputeDamageBreakageStressCycleSim
        alpha_in = alpha_in
        B_in = B_in
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        initial_alpha = initial_alpha
        # output_properties = 'eps_p eps_e eps_total I1 sts_total'
        outputs = exodus
        block = 0
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [nonADdensity]
        type = GenericConstantMaterial
        prop_names = 'nonADdensity'
        prop_values = '2700'
    []
    [./static_initial_stress_tensor]
        type = ADGenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_stress_xx        func_stress_xy     func_stress_xz 
                            func_stress_xy        func_stress_yy     func_stress_yz
                            func_stress_xz        func_stress_yz     func_stress_zz'
    [../]
    [./static_initial_strain_tensor]
        type = ADGenericFunctionRankTwoTensor
        tensor_name = static_initial_strain_tensor
        tensor_functions = 'func_strain_xx        func_strain_xy     func_strain_xz 
                            func_strain_xy        func_strain_yy     func_strain_yz
                            func_strain_xz        func_strain_yz     func_strain_zz'
    [../]
    #non-damage elastic block (id = 1)
    [./nondamagestress]
        type = ADComputeLinearElasticStress
        block = 1
    [../]
    [./elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = 30e9
        shear_modulus = 30e9
        block = 1
    [../]         
[]  

[Functions]
    #func_stress
    [func_stress_xx]
        type = ConstantFunction
        value = -50e6
    [../]
    [func_stress_xy]
        type = ConstantFunction
        value = 0
        # type = InitialStressAD
    [../]
    # [func_stress_xy]
    #     type = ParsedFunction
    #     expression = '1e4 * t'
    # []
    [func_stress_yy]
        type = ConstantFunction
        value = -50e6
        # type = InitialNormalStressAD
    [../]
    [func_stress_xz]
        type = ConstantFunction
        value = 0
    [../]
    [func_stress_yz]
        type = ConstantFunction
        value = 0
    [../]
    [func_stress_zz]
        type = ConstantFunction
        value = 0
    [../]
    #func_strain
    [func_strain_xx]
        type = ConstantFunction
        value = 0
    [../]
    [func_strain_xy]
        type = ConstantFunction
        value = 0
    [../]
    [func_strain_yy]
        type = ConstantFunction
        value = 0
    [../]
    [func_strain_xz]
        type = ConstantFunction
        value = 0
    [../]
    [func_strain_yz]
        type = ConstantFunction
        value = 0
    [../]
    [func_strain_zz]
        type = ConstantFunction
        value = 0
    [../]
    [func_initial_alpha]
        type = InitialAlphaAD
    []
    [func_top_bc]
        type = ParsedFunction
        expression = '1e-9*t'
    []
    [func_bot_bc]
        type = ParsedFunction
        expression = '-1e-9*t'
    []
    #function for adaptivity
    [solution]
        type = ConstantFunction
        value = 0
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
    solve_type = 'PJFNK'
    start_time = 0
    end_time = 1e10
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 10
    nl_abs_tol = 1e-8
    timestep_tolerance = 1e-6
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1
        cutback_factor_at_failure = 0.1
        growth_factor = 2
        enable = true
    []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]  

[Outputs]
    exodus = true
    interval = 10
    # show = 'alpha_in B_in xi_old mu_old disp_x disp_y vel_x vel_y'    
[]

[BCs]
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = right
    []
    #
    [bc_load_top_x]
        type = PresetDisplacement
        variable = disp_x
        function = func_top_bc
        boundary = top
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
    []
    [bc_fix_top_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = top
    []
    [bc_load_bottom_x]
        type = PresetDisplacement
        variable = disp_x
        function = func_bot_bc
        boundary = bottom
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
    []
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'cyclesim_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
        sub_cycling = true
    [../]
[]
  
[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub'
        variable = 'alpha_in B_in alpha_grad_x alpha_grad_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
    #we actually don't need to pass alpha and B
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'alpha_damagedvar_out B_out alpha_damagedvar_out B_out xi_old I2_old mu_old lambda_old gamma_old'
        variable = 'alpha_old B_old alpha_sub B_sub xi_old I2_old mu_old lambda_old gamma_old'
        execute_on = 'TIMESTEP_BEGIN'
    []
    #Borrowed from LevelSet module
    # [marker_to_sub]
    #     type = LevelSetMeshRefinementTransfer
    #     to_multi_app = sub_app
    #     source_variable = thresholdmarker
    #     variable = thresholdmarker
    #     check_multiapp_execute_on = false
    # []
[]

# [Adaptivity]
#     marker = thresholdmarker
#     max_h_level = 2
#     [Markers]
#         [thresholdmarker]
#             type = ValueThresholdMarker
#             refine = 0.5
#             variable = B_in
#         []
#     []
# []

# #this is added only to activate "LevelSetMeshRefinementTransfer"
# [Problem]
#     type = LevelSetProblem
# []