[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = './meshfile/tpv2052dm.msh'
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
[]

[GlobalParams]

    displacements = 'disp_x disp_y'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 10e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 10e9
    
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
    Cd_constant = 300

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 500

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.01 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.5
    
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
    [xi_computed]
        order = CONSTANT
        family = MONOMIAL
    []
    [initial_damage_aux]
        order = CONSTANT
        family = MONOMIAL     
    []
    #
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
[]

[AuxKernels]
    [compute_xi]
        type = CompXi3D
        variable = xi_computed
        execute_on = 'TIMESTEP_END'
    []
    #
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
    # [get_initial_damage_aux]   
    #     type = SolutionAux
    #     variable = initial_damage_aux
    #     solution = init_sol_components
    #     from_variable = initial_damage
    # []
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
        use_displaced_mesh = false
        variable = disp_x
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2700'
    []   
    [nonADdensity]
        type = GenericConstantMaterial
        prop_names = 'nonADdensity'
        prop_values = '2700'
    []  
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar'
        outputs = exodus
        block = 1
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain deviatroic_stress'
        outputs = exodus
        block = 1
    []
    # elastic
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 1e10
        shear_modulus = 1e10
        block = 2
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
        block = 2
    []
    # [initialdamage]
    #     type = ParsedMaterial
    #     property_name = initial_damage
    #     coupled_variables = initial_damage_aux
    #     expression = 'initial_damage_aux'
    #     outputs = exodus
    # []
    [initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = 0
    []  
[]  

[Functions]
    [func_top_bc]
        type = ParsedFunction
        expression = '1e-9*t'
    []
    [func_bot_bc]
        type = ParsedFunction
        expression = '-1e-9*t'
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
    start_time = 0
    end_time = 1e100
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 20
    nl_abs_tol = 1e-8
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    # dt = 20
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 0.01
        cutback_factor_at_failure = 0.5
        optimal_iterations = 10
        growth_factor = 1.5
        enable = true
        reject_large_step_threshold = 0.01
        reject_large_step = true
    []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

[Outputs] 
    exodus = true
    time_step_interval = 100
    [./my_checkpoint]
        type = Checkpoint
        num_files = 2
        interval = 100
    [../]
[]

[BCs]
    [bc_load_top_x]
        type = FunctionDirichletBC
        variable = disp_x
        function = func_top_bc
        boundary = top
    []
    [bc_fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = bottom
    []
    [bc_fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    [bc_fix_left_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = left
    []
    [bc_fix_right_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = right
    []
    [./neumann_top_y]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -120e6
    [../]
    [./neumann_left_x]
        type = NeumannBC
        variable = disp_x
        boundary = left
        value = 135e6
    [../]
    [./neumann_right_x]
        type = NeumannBC
        variable = disp_x
        boundary = right
        value = -135e6
    [../]
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = './static_solve_out.e'
      system_variables = 'disp_x disp_y'
      timestep = LATEST
      force_preaux = true
    [../]
[]

# [ICs]
#     [disp_x_ic]
#       type = SolutionIC
#       variable = disp_x
#       solution_uo = init_sol_components
#       from_variable = disp_x
#     []
#     [disp_y_ic]
#       type = SolutionIC
#       variable = disp_y
#       solution_uo = init_sol_components
#       from_variable = disp_y
#     []
# []

[Problem]
    restart_file_base = dynamic_solve_my_checkpoint_cp/LATEST  # You may also use a specific number here
[]