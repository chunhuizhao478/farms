[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/cdbm_tpv2053d_buried.msh'
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
        coord = '-20000  -20000  -20000;
                  20000  -20000  -20000;
                 -20000  -20000   20000;
                  20000  -20000   20000'
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
    Cd_constant = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

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
    [initial_damage_aux]
        order = CONSTANT
        family = MONOMIAL
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
    [initial_damage]
        type = SolutionAux
        variable = initial_damage_aux
        solution = init_sol_components
        from_variable = initial_damage
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

[AuxKernels]
[]

[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        # outputs = exodus
    [] 
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2700
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3Ddebug
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi'
        outputs = exodus
    [] 
    [initial_damage]
        type = ParsedMaterial
        property_name = initial_damage
        coupled_variables = initial_damage_aux
        expression = 'initial_damage_aux'
        outputs = exodus
    []
    # [damage_perturb]
    #     type = DamagePerturbationSperical
    #     nucl_center = '0 -2500 0'
    #     e_damage = 0.3
    #     e_sigma = 1e3
    #     duration = 1e-1
    #     outputs = exodus
    # []
    [damage_perturb]
        type = DamagePerturbationSquare
        nucl_center = '0 -10000 0'
        e_damage = 0.3
        thickness = 100
        length = 400
        duration = 1e-1
        outputs = exodus
    []
[]  

[Functions]
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve_buried/static_solve_out.e'
      system_variables = 'disp_x disp_y disp_z initial_damage'
      timestep = LATEST
      force_preaux = true
    [../]
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]
  
[Executioner]
    type = Transient
    dt = 1e-4
    end_time = 50.0
    # num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]

[Outputs] 
    exodus = true
    time_step_interval = 1000
    [sample_snapshots]
        type = Exodus
        time_step_interval = 2000
    []
    [./checkpoint]
        type = Checkpoint
        wall_time_interval = 4000 # interval length in seconds
    [../]    
[]

#We assume the simulation is loaded with compressive pressure and shear stress
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
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        value = 55e6
    []
    [pressure_shear_back]
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        value = -55e6   
    []
    [pressure_shear_left]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        value = -55e6
    []
    [pressure_shear_right]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        value = 55e6     
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