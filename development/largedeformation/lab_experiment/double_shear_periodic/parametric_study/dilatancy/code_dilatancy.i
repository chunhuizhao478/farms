[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 100
        ny = 20
        xmin = 0
        xmax = 0.05
        ymin = 0
        ymax = 0.01
    [] 
    [./box]
        type = SubdomainBoundingBoxGenerator
        input = msh
        block_id = 1
        bottom_left = '0 0 0'
        top_right = '0.05 0.004 0'
    []
    [./box2]
        type = SubdomainBoundingBoxGenerator
        input = box
        block_id = 0
        bottom_left = '0 0.004 0'
        top_right = '0.05 0.006 0'
    []
    [./box3]
        type = SubdomainBoundingBoxGenerator
        input = box2
        block_id = 2
        bottom_left = '0 0.006 0'
        top_right = '0.05 0.01 0'
    []
    [internal_top]
        type = SideSetsBetweenSubdomainsGenerator
        new_boundary = 'internal_top'
        primary_block = 1
        paired_block = 0
        input = box3
    []  
    [internal_bottom]
        type = SideSetsBetweenSubdomainsGenerator
        new_boundary = 'internal_bottom'
        primary_block = 2
        paired_block = 0
        input = internal_top
    []  
[]

[GlobalParams]

    displacements = 'disp_x disp_y'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 32e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 32e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.75
    
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
    Cd_constant = 1.86e5

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.05 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-11 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
    
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
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
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
        outputs = exodus
        block = 0
        #add dilatancy
        add_dilatancy_compaction_anand = true
        anand_param_go = 0.04
        anand_param_eta_cv = 0.1
        anand_param_p = 0.8
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
        lambda = 64e9
        shear_modulus = 64e9
        block = '1 2'
    []
    [compute_stress]
        type = ComputeStVenantKirchhoffStress
        large_kinematics = true
        block = '1 2'
    []
[]  

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '1e-6 * t'
    []
    [applied_strain_top]
        type = ParsedFunction
        expression = '1e-6 * t / 0.05'
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
    end_time = 150 #extend the time
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 5
    nl_abs_tol = 1e-8
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    # dt = 0.1
    [./TimeIntegrator]
        type = ImplicitEuler
        # type = BDF2
        # type = CrankNicolson
    [../]
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-3
        cutback_factor_at_failure = 0.5
        optimal_iterations = 8
        growth_factor = 1.1
        max_time_step_bound = 0.1
    []
[]

[Outputs] 
    exodus = true
    time_step_interval = 1
    [./csv]
        type = CSV
        time_step_interval = 1
        show = 'strain_x react_x'
    [../]
[]

[BCs]
    #fix bottom
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
    #load on top
    [load_top_y2]
        type = NeumannBC
        variable = disp_y
        boundary = top
        value = -50e6
    []
    #periodic boundary
    [./Periodic]
        [./x]
          variable = disp_x
          primary = left
          secondary = right
          translation = '0.05 0 0'
        [../]
        [./y]
          variable = disp_y
          primary = left
          secondary = right
          translation = '0.05 0 0'
        [../]
    [../]
    #displacement rate    
    [applied_top_x2]
        type = FunctionDirichletBC
        variable = disp_x
        boundary = top
        function = applied_load_top
    [] 
    #fix lateral boundaries
    [fix_left_y]
        type = DirichletBC
        variable = disp_y
        boundary = left
        value = 0
    [] 
    [fix_right_y]
        type = DirichletBC
        variable = disp_y
        boundary = right
        value = 0
    []
    #fix internal boundaries
    [fix_internal_top_y]
        type = DirichletBC
        variable = disp_y
        boundary = internal_top
        value = 0
    []
    [fix_internal_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = internal_bottom
        value = 0
    []
[]

#compute the reaction force on the top boundary
[Postprocessors]
    [./react_x]
      type = SidesetReaction
      direction = '1 0 0'
      stress_tensor = pk2_stress
      boundary = top
    [../]
    [./strain_x]
        type = FunctionValuePostprocessor
        function = applied_strain_top
    []
[]

