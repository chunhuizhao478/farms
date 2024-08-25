[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/cdbm_tpv2053d_small.msh'
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
        new_boundary = 'left right bottom top back front'
    []    
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '-8000  -8000  -8000;
                  8000  -8000  -8000;
                 -8000  -8000   8000;
                  8000  -8000   8000'
        new_boundary = corner_ptr
        input = sidesets
    [] 
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
    xi_d = -0.9
  
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 1e5

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-8
  
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

    #
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [vel_x]
    []  
    [vel_y]
    []
    [vel_z]
    []
[]

[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_z
      coupled = disp_z
      execute_on = 'TIMESTEP_END'
    []
[]
  
[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
    [./inertia_x]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y
    []
    [./inertia_z]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_z
    []
[]
  
[Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeDamageBreakageStress3D
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi epsilon_eq'
        outputs = exodus
    []
    [strain]
        type = ComputeSmallStrain
        outputs = exodus
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2700
    []
    [nonADdensity]
        type = GenericConstantMaterial
        prop_names = nonADdensity
        prop_values = 2700
    []
    [./static_initial_stress_tensor]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_stress_xx     func_stress_xy      func_stress_xz 
                            func_stress_xy     func_stress_yy      func_stress_yz
                            func_stress_xz     func_stress_yz      func_stress_zz'
    [../]
    [./static_initial_strain_tensor]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_strain_tensor
        tensor_functions = 'func_strain_xx     func_strain_xy      func_strain_xz 
                            func_strain_xy     func_strain_yy      func_strain_yz
                            func_strain_xz     func_strain_yz      func_strain_zz'
    [../]
    [./I1_initial]
        type = GenericFunctionMaterial
        prop_names = I1_initial
        prop_values = func_I1  
    []
    [./I2_initial]
        type = GenericFunctionMaterial
        prop_names = I2_initial
        prop_values = func_I2 
    []  
    [./xi_initial]
        type = GenericFunctionMaterial
        prop_names = xi_initial
        prop_values = func_xi
    []  
    [./initial_damage]
        type = GenericFunctionMaterial
        prop_names = initial_damage
        prop_values = func_initial_damage
    [] 
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve_small/static_solve_out.e'
      system_variables = 'disp_x disp_y disp_z I1_initial I2_initial xi_initial initial_damage mechanical_strain_00 mechanical_strain_01 mechanical_strain_02 mechanical_strain_11 mechanical_strain_12 mechanical_strain_22 stress_00 stress_01 stress_02 stress_11 stress_12 stress_22'
      timestep = LATEST
      force_preaux = true
    [../]
[]
  
[Functions]
    [func_strain_xx]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_00
      execute_on = 'INITIAL'
    [../]
    [func_strain_xy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_01
      execute_on = 'INITIAL'
    [../]
    [func_strain_xz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_02
      execute_on = 'INITIAL'
    [../]
    [func_strain_yy]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_11
      execute_on = 'INITIAL'
    [../]
    [func_strain_yz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_12
      execute_on = 'INITIAL'
    [../]
    [func_strain_zz]
      type = SolutionFunction
      solution = init_sol_components
      from_variable = mechanical_strain_22
      execute_on = 'INITIAL'
    []
    #
    [func_stress_xx]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_00
        execute_on = 'INITIAL'
    [../]
    [func_stress_xy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_01
        execute_on = 'INITIAL'
    [../]
    [func_stress_xz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_02
        execute_on = 'INITIAL'
    [../]
    [func_stress_yy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_11
        execute_on = 'INITIAL'
    [../]
    [func_stress_yz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_12
        execute_on = 'INITIAL'
    [../]
    [func_stress_zz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = stress_22
        execute_on = 'INITIAL'
    [../] 
    #
    [func_I1]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = I1_initial
        execute_on = 'INITIAL'
    [../] 
    [func_I2]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = I2_initial
        execute_on = 'INITIAL'
    [../] 
    [func_xi]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = xi_initial
        execute_on = 'INITIAL'
    [../]   
    [func_initial_damage]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = initial_damage
        execute_on = 'INITIAL'
    [../]   
[]

#0.4/5773
[Executioner]
    type = Transient
    dt = 1e-4
    end_time = 1.0
    # num_steps = 8000
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]
  
[Outputs]
    exodus = true
    time_step_interval = 100
    [sample_snapshots]
        type = Exodus
        interval = 200
    []
    [snapshots]
        type = Exodus
        interval = 50
        overwrite = true
    []
[]

[BCs]
    [pressure_right]
        type = Pressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        factor = 135e6
    []
    [pressure_left]
        type = Pressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        factor = 135e6
    []
    [pressure_front]
        type = Pressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        factor = 120e6
    []
    [pressure_back]
        type = Pressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        factor = 120e6        
    []
    [pressure_top]
        type = Pressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = top
        factor = 127.5e6         
    []
    [pressure_bottom]
        type = Pressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = bottom
        factor = 127.5e6              
    []
    #
    [pressure_shear_front]
        type = NeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        value = 65e6
    []
    [pressure_shear_back]
        type = NeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        value = -65e6   
    []
    [pressure_shear_left]
        type = NeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        value = -65e6
    []
    [pressure_shear_right]
        type = NeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        value = 65e6     
    []
    #
    [fix_ptr_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_z]
        type = DirichletBC
        variable = disp_z
        value = 0
        boundary = corner_ptr
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
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]