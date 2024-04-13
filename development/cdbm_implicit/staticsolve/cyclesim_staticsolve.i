#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/cyclesim.msh'
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
    C_g = 1e-5
  
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
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
    []
[]

[AuxKernels]
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
    [get_alpha_old]
        type = ADMaterialRealAux
        property = alpha_damagedvar
        variable = alpha_damagedvar_out
        execute_on = 'INITIAL TIMESTEP_BEGIN'   
    []
    [get_B_old]
        type = ADMaterialRealAux
        property = B
        variable = B_out
        execute_on = 'INITIAL TIMESTEP_BEGIN'   
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
        function = func_stress_xy
        variable = shear_stress_applied
    []
[]

[Materials]
    [damagestress]
        type = ADComputeDamageBreakageStressDev
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
    [func_stress_yy]
        type = ConstantFunction
        value = -50e6
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
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_factor_shift_type'
      petsc_options_value = 'lu       NONZERO'
    []
[]
  
[Executioner]
    type = Steady
    solve_type = Newton
    l_max_its = 10
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-10
    automatic_scaling = true
    line_search = 'none'
[]  

[Outputs]
    exodus = true
[]

[BCs]
    #pressure
    ##compression
    [pressure_top]
        type = ADPressure
        variable = disp_y
        boundary = top
        factor = 50e6
    []
    [pressure_bottom]
        type = ADPressure
        variable = disp_y
        boundary = bottom
        factor = 50e6
    []
    [pressure_left]
        type = ADPressure
        variable = disp_x
        boundary = left
        factor = 50e6
    []
    [pressure_right]
        type = ADPressure
        variable = disp_x
        boundary = right
        factor = 50e6
    []  
[]