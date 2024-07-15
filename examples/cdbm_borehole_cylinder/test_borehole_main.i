#Borehole Breakouts
#Sierra white granite
#density: 2640 kg/m^3
#youngs modulus: 67 GPa
#internal friction angle: 46
#poisson's ratio: 0.32

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  './meshfile/cylinder_w_hole.msh'
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
  
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 4.083e9
      
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 6.125e9
  
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.968
  
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.968
  
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
    gamma_damaged_r = 6.728e9
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.746
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.75
    a0 = 1.166e9
    a1 = -3.495e9
    a2 = 3.082e9
    a3 = -0.660e9
  
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
        order = FIRST
        family = LAGRANGE
    []
    [./B_in]
        order = FIRST
        family = LAGRANGE
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
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
    []
    [dispkernel_z]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        use_displaced_mesh = true
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_z]
        type = ADInertialForce
        variable = disp_z
        velocity = vel_z
        acceleration = accel_z
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
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
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = timestep_end
    []
    #
    [get_xi_old]
        type = ADMaterialRealAux
        property = xi
        variable = xi_old
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_I2_old]
        type = ADMaterialRealAux
        property = I2
        variable = I2_old
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_mu_old]
        type = ADMaterialRealAux
        property = shear_modulus
        variable = mu_old
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_lambda_old]
        type = ADMaterialRealAux
        property = lambda
        variable = lambda_old
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_gamma_old]
        type = ADMaterialRealAux
        property = gamma_damaged
        variable = gamma_old
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_alpha_old]
        type = ADMaterialRealAux
        property = alpha_damagedvar
        variable = alpha_damagedvar_out
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [get_B_old]
        type = ADMaterialRealAux
        property = B
        variable = B_out
        execute_on = 'INITIAL TIMESTEP_END'
    []
    #principal strain
    [get_principal_strain]
        type = ADMaterialRateRealAux
        property = principal_strain
        variable = principal_strain_rate_out
        execute_on = 'INITIAL TIMESTEP_END'
    []
    [record_applied_shear_stress]
        type = FunctionAux
        function = func_stress_xz
        variable = shear_stress_applied
        execute_on = 'INITIAL'
    []
[]

[Materials]
    [damagestress]
        type = ADComputeDamageBreakageStress3D
        alpha_in = alpha_in
        B_in = B_in
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        initial_alpha = initial_alpha
        # output_properties = 'eps_p eps_e eps_total I1 sts_total'
        outputs = exodus
        # outputs = nemesis
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '1918.77'
    []
    [nonADdensity]
        type = GenericConstantMaterial
        prop_names = 'nonADdensity'
        prop_values = '1918.77'
    []
    #we assume symmetry here
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
[]  

[Functions]
    #func_stress
    #left-right
    [func_stress_xx]
        type = ConstantFunction
        value = -70e6
    [../]
    [func_stress_xy]
        type = ConstantFunction
        value = 0
    [../]
    #
    [func_stress_yy]
        type = ConstantFunction
        value = -20e6
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
        value = -30e6
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
    #
    [func_bc]
        type = ParsedFunction
        expression = '1e6 * t'
    []
[]

[BCs]
    #
    [disp_top]
        type = Pressure
        boundary = top
        variable = disp_y
        factor = 20e6
    []
    [disp_bottom]
        type = Pressure
        boundary = bottom
        variable = disp_y
        factor = 20e6
    []
    #
    [disp_left]
        type = Pressure
        boundary = left
        variable = disp_x
        factor = 70e6
    []
    [disp_right]
        type = Pressure
        boundary = right
        variable = disp_x
        factor = 70e6
    []
    #
    [disp_front]
        type = Pressure
        boundary = front
        variable = disp_z
        factor = 30e6
    []
    [disp_back]
        type = Pressure
        boundary = back
        variable = disp_z
        factor = 30e6
    []
    #
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr1
        value = 0
    []
    [./fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr1
        value = 0
    []
    [./fix_cptr1_z]
        type = DirichletBC
        variable = disp_z
        boundary = corner_ptr1
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
    solve_type = 'NEWTON'
    start_time = 0
    end_time = 18000
    # num_steps = 10
    l_max_its = 20
    nl_rel_tol = 1e-6
    nl_max_its = 10
    nl_abs_tol = 1e-8
    timestep_tolerance = 1e-6
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    # dt = 1
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1e-3
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
    interval = 1
    show = 'alpha_in B_in xi_old mu_old disp_x disp_y disp_z vel_x vel_y vel_z'    
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'test_borehole_sub.i'
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
        source_variable = 'alpha_in B_in alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old'
        variable = 'alpha_old B_old alpha_sub B_sub xi_old I2_old mu_old lambda_old gamma_old'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]