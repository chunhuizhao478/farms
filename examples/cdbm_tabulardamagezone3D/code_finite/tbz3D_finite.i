#implicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/tabulardamagezone_small.msh'
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
        coord = '-100 -100 -100'
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
    xi_d = -0.9
  
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 10

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
    C_g = 1e-15
  
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
[]

[Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
        static_initialization = true
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
        static_initialization = true
    []
    [dispkernel_z]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        use_displaced_mesh = true
        static_initialization = true
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
[]

[Materials]
    [damagestress]
        type = ADComputeFiniteStrainDamageBreakageStress
        extra_stress_names = stress_function
        outputs = exodus
    []
    [strain]
        type = ADComputeFiniteStrain
        displacements = 'disp_x disp_y disp_z'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []
    [nonADdensity]
        type = ADGenericConstantMaterial
        prop_names = 'nonADdensity'
        prop_values = '2700'
    []
    [./strain_from_initial_strain]
        type = ADComputeEigenstrainFromSolution
        initial_strain = 'func_strain_xx func_strain_xy func_strain_xz  
                          func_strain_xy func_strain_yy func_strain_yz 
                          func_strain_xz func_strain_yz func_strain_zz'
        eigenstrain_name = ini_strain
    [../]
    [./stress_function]
        type = GenericFunctionRankTwoTensor
        tensor_name = stress_function
        tensor_functions = 'func_stress_xx func_stress_xy func_stress_xz  
                            func_stress_xy func_stress_yy func_stress_yz 
                            func_stress_xz func_stress_yz func_stress_zz'
    [../]
    [initialdamage]
        type = ADInitialDamage
    []     
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = 'disp_x disp_y disp_z'
      timestep = LATEST
      force_preaux = true
    [../]
    #add initial strain components
    [load_strain_xx]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_00
      timestep = LATEST
      force_preaux = true
    []
    [load_strain_xy]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_01
      timestep = LATEST
      force_preaux = true
    []
    [load_strain_xz]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_02
      timestep = LATEST
      force_preaux = true
    []
    [load_strain_yy]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_11
      timestep = LATEST
      force_preaux = true
    [] 
    [load_strain_yz]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_12
      timestep = LATEST
      force_preaux = true
    []  
    [load_strain_zz]
      type = SolutionUserObject
      mesh = ../static_solve/static_solve_out.e
      system_variables = mechanical_strain_22
      timestep = LATEST
      force_preaux = true
    []  
    #add initial stress components
    [load_stress_xx]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_00
        timestep = LATEST
        force_preaux = true
    []
    [load_stress_xy]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_01
        timestep = LATEST
        force_preaux = true
    []
    [load_stress_xz]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_02
        timestep = LATEST
        force_preaux = true
    []
    [load_stress_yy]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_11
        timestep = LATEST
        force_preaux = true
    [] 
    [load_stress_yz]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_12
        timestep = LATEST
        force_preaux = true
    []  
    [load_stress_zz]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = stress_22
        timestep = LATEST
        force_preaux = true
    [] 
[]

[Functions]
     [func_strain_xx]
      type = SolutionFunction
      solution = load_strain_xx
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_strain_xy]
      type = SolutionFunction
      solution = load_strain_xy
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_strain_xz]
      type = SolutionFunction
      solution = load_strain_xz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_strain_yy]
      type = SolutionFunction
      solution = load_strain_yy
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_strain_yz]
      type = SolutionFunction
      solution = load_strain_yz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_strain_zz]
      type = SolutionFunction
      solution = load_strain_zz
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../] 
    #
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
    [func_stress_xz]
        type = SolutionFunction
        solution = load_stress_xz
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_stress_yy]
        type = SolutionFunction
        solution = load_stress_yy
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_stress_yz]
        type = SolutionFunction
        solution = load_stress_yz
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
    [func_stress_zz]
        type = SolutionFunction
        solution = load_stress_zz
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]   
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
      petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '      
    []
  []
  
[Executioner]
    type = Transient
    solve_type = PJFNK
    start_time = 0
    end_time = 800
    num_steps = 1000
    l_max_its = 10
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 10
    nl_abs_tol = 1e-8
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1e-4
    [TimeIntegrator]
        type = NewmarkBeta
    []
[]  

[Outputs]
    exodus = true
    interval = 10
    # show = 'alpha B disp_x disp_y disp_z vel_x vel_y vel_z'
    print_linear_residuals=true
[]

[BCs]
    [./dashpot_top_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = top
    []
    [./dashpot_top_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = top
    []
    [./dashpot_top_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = top
    []
    [./dashpot_bottom_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = bottom
    []
    [./dashpot_bottom_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = bottom
    []
    [./dashpot_left_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = left
    []
    [./dashpot_left_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = left
    []
    [./dashpot_left_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = left
    []
    [./dashpot_right_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = right
    []
    [./dashpot_right_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = right
    []
    [./dashpot_right_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = right
    []
    [./dashpot_back_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = back
    []
    [./dashpot_back_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = back
    []
    [./dashpot_back_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = back
    []
    [./dashpot_front_x]
        type = ADNonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = front
    []
    [./dashpot_front_y]
        type = ADNonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = front
    []
    [./dashpot_front_z]
        type = ADNonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = 5773.50
        shear_wave_speed = 3333.33
        boundary = front
    []
    #
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