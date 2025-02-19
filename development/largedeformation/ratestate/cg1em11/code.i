[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 3
        nx = 10
        ny = 10
        nz = 10
        xmin = 0
        xmax = 1
        ymin = 0
        ymax = 1
        zmin = 0
        zmax = 1
    [] 
    [./box]
        type = SubdomainBoundingBoxGenerator
        input = msh
        block_id = 1
        bottom_left = '0 0 0'
        top_right = '1 0.5 1'
    []
    [./box2]
        type = SubdomainBoundingBoxGenerator
        input = box
        block_id = 1
        bottom_left = '0 0.6 0'
        top_right = '1.0 1.0 1.0'
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
    Cd_constant = 30

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 500

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.01 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-11 #
    
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
    [disp_z]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    [xi_computed]
        order = CONSTANT
        family = MONOMIAL
    []
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    [compute_xi]
        type = CompXi3D
        variable = xi_computed
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    #
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    []
    [Vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
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
[]

[Materials]
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        output_properties = 'deformation_gradient'
        outputs = exodus
    []
    # damage
    [damage_mat]
        type = DamageBreakageMaterial
        output_properties = 'alpha_damagedvar B_damagedvar velgrad_L'
        # use_vels_build_L = true
        # vel_x = vel_x
        # vel_y = vel_y
        # vel_z = vel_z
        outputs = exodus
        block = 0
    [] 
    [initial_damage]
        type = GenericConstantMaterial
        prop_names = 'initial_damage'
        prop_values = 0
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Debug
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain plastic_deformation_gradient_det'
        outputs = exodus
        block = 0
    []
    # elastic
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 30e9
        shear_modulus = 30e9
        block = 1
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        block = 1
    []
[]  

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '1e-4 * t'
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
    end_time = 400
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 5
    nl_abs_tol = 1e-6
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    dt = 1
    [./TimeIntegrator]
        # type = ImplicitEuler
        type = BDF2
        # type = CrankNicolson
    [../]
[]

[Outputs] 
    exodus = true
    time_step_interval = 1
[]

[BCs]
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = bottom
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = bottom
        value = 0
    []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []   
    [fix_top_y]
        type = DirichletBC
        variable = disp_y
        boundary = top
        value = 0
    [] 
    [applied_top_x]
        type = FunctionDirichletBC
        variable = disp_x
        boundary = top
        function = applied_load_top
    []
    [./Pressure]
        [static_pressure_back]
            boundary = back
            factor = 50e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []  
        [static_pressure_front]
            boundary = front
            factor = 50e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []    
        [static_pressure_left]
            boundary = left
            factor = 50e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []  
        [static_pressure_right]
            boundary = right
            factor = 50e6
            displacements = 'disp_x disp_y disp_z'
            use_displaced_mesh = false
        []         
    []    
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = 'initial_load_check_out.e'
      system_variables = 'disp_x disp_y disp_z'
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
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]