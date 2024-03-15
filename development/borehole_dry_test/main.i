#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  './meshfile/borehole_wofaults.msh'
    []
[]

[GlobalParams]
  
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'

    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 14e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 14e9

    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.929

    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.929

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
    gamma_damaged_r = 15.549e9

    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.7631

    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m

    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#

    # #coefficients
    # chi = 0.75
    a0 = 2.788e9
    a1 = -8.6937e9
    a2 = 8.4215e9
    a3 = -2.3806e9

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
[]

[Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
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
    []
    [get_I2_old]
        type = ADMaterialRealAux
        property = I2
        variable = I2_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_mu_old]
        type = ADMaterialRealAux
        property = shear_modulus
        variable = mu_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_lambda_old]
        type = ADMaterialRealAux
        property = lambda
        variable = lambda_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_gamma_old]
        type = ADMaterialRealAux
        property = gamma_damaged
        variable = gamma_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    # #initial alpha
    # [initial_alpha]
    #     type = FunctionAux
    #     variable = alpha_in
    #     function = func_initial_alpha
    #     execute_on = 'INITIAL'
    # []
[]

[Materials]
    # [Elasticity_tensor]
    #     type = ADComputeElasticityTensor
    #     fill_method = symmetric_isotropic
    #     C_ijkl = '210e9 0'
    # []
    # [stress]
    #     type = ADComputeLinearElasticStress
    # []
    [damagestress]
        type = ADComputeDamageBreakageStressv2
        alpha_in = alpha_in
        B_in = B_in
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        # output_properties = 'eps_p eps_e eps_total I1 sts_total'
        outputs = exodus
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2650'
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
[]  

[Functions]
    #func_stress
    [func_stress_xx]
        type = SolutionFunction
        solution = load_stress_xx
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_stress_xy]
        type = SolutionFunction
        solution = load_stress_xy
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_stress_yy]
        type = SolutionFunction
        solution = load_stress_yy
        execute_on = 'INITIAL TIMESTEP_BEGIN'
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
    # [func_initial_alpha]
    #     type = InitialAlphaAD
    # []
[]

[UserObjects]
    [load_stress_xx]
      type = SolutionUserObject
      mesh = ./static_solve/borehole_static_solve_out.e
      system_variables = stress_xx_saved
      execute_on = 'INITIAL TIMESTEP_BEGIN'
      timestep = LATEST
      force_preaux = true
    []
    [load_stress_xy]
      type = SolutionUserObject
      mesh = ./static_solve/borehole_static_solve_out.e
      system_variables = stress_xy_saved
      execute_on = 'INITIAL TIMESTEP_BEGIN'
      timestep = LATEST
      force_preaux = true
    []
    [load_stress_yy]
      type = SolutionUserObject
      mesh = ./static_solve/borehole_static_solve_out.e
      system_variables = stress_yy_saved
      execute_on = 'INITIAL TIMESTEP_BEGIN'
      timestep = LATEST
      force_preaux = true
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
    end_time = 800
    num_steps = 200
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-4
    nl_max_its = 10
    nl_abs_tol = 1e-6
    timestep_tolerance = 1e-6
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1
[]  

[Outputs]
    exodus = true
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'sub.i'
        execute_on = 'TIMESTEP_BEGIN'
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
        source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old'
        variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]