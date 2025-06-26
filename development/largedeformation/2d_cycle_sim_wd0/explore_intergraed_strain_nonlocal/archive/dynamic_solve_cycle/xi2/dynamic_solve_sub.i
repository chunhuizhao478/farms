#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../../mesh/mesh_local.msh'
    []
    [mergeblock1]
        type = RenameBlockGenerator
        input = msh
        old_block = '1 3'
        new_block = '4 4'
    []
    [leftstrength]
        type = SubdomainBoundingBoxGenerator
        block_id = 5
        bottom_left = '-6000 -2000 0'
        top_right = '-2000 2000 0'
        input = mergeblock1
    []
    [rightstrength]
        type = SubdomainBoundingBoxGenerator
        block_id = 5
        bottom_left = '2000 -2000 0'
        top_right = '6000 2000 0'
        input = leftstrength
    []
    [./sidesets]
        input = rightstrength
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '-30000 -30000 0'
        new_boundary = corner_ptr
        input = sidesets
    []
    [./extranodeset2]
        type = ExtraNodesetGenerator
        coord = '30000 -30000 0'
        new_boundary = corner_ptr
        input = extranodeset1
    []
[]

[GlobalParams]
    
    ##----continuum damage breakage model----##
    #initial lambda value (FIRST lame constant) [Pa]
    lambda_o = 32.04e9
        
    #initial shear modulus value (FIRST lame constant) [Pa]
    shear_modulus_o = 32.04e9
    
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

    #if option 2, use Cd_constant #specify by auxiliary variable
    Cd_constant = -1.0

    #strain rate dependent Cd options
    m_exponent = 0.8
    strain_rate_hat = 1e-8
    cd_hat = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd #specify by auxiliary variable
    CdCb_multiplier = 100

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf" #specify by auxiliary variable
    C_1 = 1e-4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3

    #diffusion parameter #close the gradient
    D_diffusion = 0

[]

#this sub-app solves damage/breakage evolution equations
[Variables]
    [alpha_damagedvar_sub]
        order = FIRST
        family = LAGRANGE     
    []
    [B_damagedvar_sub]
        order = FIRST
        family = LAGRANGE    
    []
[]

[AuxVariables]
    [xi_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    [I2_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_damage_sub_aux]
        order = FIRST
        family = LAGRANGE
    []
    #deviatroic_strain_rate
    [deviatroic_strain_rate_sub_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [bounds_dummy]
        order = FIRST
        family = LAGRANGE
    []
    #
    [structural_stress_coefficient_sub]
        order = FIRST
        family = MONOMIAL
    []
    #
    [Cd_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xio_aux]
        order = FIRST
        family = LAGRANGE
    []
    [xid_aux]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    #damagevar
    [time_derivative_alpha]
        type = TimeDerivative
        variable = alpha_damagedvar_sub
    []
    [diffusion_alpha]
        type = DamageEvolutionDiffusion
        variable = alpha_damagedvar_sub
        coupled = B_damagedvar_sub
        block = '4 5'
    []
    [forcing_term_alpha]
        type = DamageEvolutionConditionalForcing
        variable = alpha_damagedvar_sub
        coupled = B_damagedvar_sub
        block = '4 5'
    []
    [perturb_source_alpha]
        type = PerturbationSource
        variable = alpha_damagedvar_sub
        damage_source = 'damage_perturbation'
        block = '4 5'
    []
    #breakagevar
    [time_derivative_B]
        type = TimeDerivative
        variable = B_damagedvar_sub
    []
    [forcing_term_B]
        type = BreakageEvolutionConditionalForcing
        variable = B_damagedvar_sub
        coupled = alpha_damagedvar_sub
        block = '4 5'
    []
[]

[Bounds]
    [alpha_damagedvar_upper_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = alpha_damagedvar_sub
      bound_type = upper
      bound_value = 1
    []
    [alpha_damagedvar_lower_bound]
      type = VariableConstantBounds
      variable = bounds_dummy
      bounded_variable = alpha_damagedvar_sub
      bound_type = lower
      bound_value = initial_damage_sub_aux
    []
    [B_damagedvar_upper_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = B_damagedvar_sub
      bound_type = upper
      bound_value = 1
    []
    [B_damagedvar_lower_bound]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = B_damagedvar_sub
      bound_type = lower
      bound_value = 0
    []
[]

[AuxKernels]
    [get_initial_damage]
        type = SolutionAux
        variable = initial_damage_sub_aux
        solution = init_sol_components
        from_variable = alpha_damagedvar_output
        execute_on = 'TIMESTEP_BEGIN'
    []
    [get_Cd]
        type = MaterialRealAux
        variable = Cd_aux
        property = Cd
    []
    #
    [get_xio]
        type = FunctionAux
        variable = xio_aux
        function = func_spatial_xio
    []
    [get_xid]
        type = FunctionAux
        variable = xid_aux
        function = func_spatial_xid
    []
[]

[Functions]
    [func_spatial_xio]
        type = SpatialDamageBreakageParameters
        W = 6e3 #half the total width
        w = 1e3
        max_val = -0.8
        min_val = 1.8
    []
    [func_spatial_xid]
        type = SpatialDamageBreakageParameters
        W = 6e3 #half the total width
        w = 1e3
        max_val = -0.9
        min_val = 1.8
    []
[]

[Materials]
    # # damage
    [damage_mat]
        type = DiffusedDamageBreakageMaterialSubApp
        I2_aux = I2_sub_aux
        xi_aux = xi_sub_aux
        initial_damage_aux = initial_damage_sub_aux
        #use strain rate dependent Cd
        use_cd_strain_dependent = true
        strain_rate = deviatroic_strain_rate_sub_aux
        #use spatial xio
        use_spatial_xio = true
        xio_aux = xio_aux
        #use spatial xid
        use_spatial_xid = true
        xid_aux = xid_aux
    []
    #add shear perturbation to the system
    [damage_perturbation]
        type = PerturbationRadialSource
        nucl_center = '0 0 0'
        peak_value = 0.3
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'damage'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
        outputs = exodus
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
    start_time = -1e-12
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 10
    nl_abs_tol = 1e-10
    petsc_options_iname = '-snes_type'
    petsc_options_value = 'vinewtonrsls'
    verbose = true
    # dt = 1e-2
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-2
        cutback_factor_at_failure = 0.5
        optimal_iterations = 8
        growth_factor = 1.1
        max_time_step_bound = 1e7
    []
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../../static_solve_newbc/static_solve_local_out.e'
      system_variables = 'alpha_damagedvar_output B_damagedvar_output'
      timestep = LATEST
      force_preaux = true
      execute_on = 'INITIAL'
    [../]
[]

[ICs]
    [alpha_damagedvar_sub_ic]
        type = SolutionIC
        variable = alpha_damagedvar_sub
        solution_uo = init_sol_components
        from_variable = alpha_damagedvar_output
    []  
    [B_damagedvar_sub_ic]
        type = SolutionIC
        variable = B_damagedvar_sub
        solution_uo = init_sol_components
        from_variable = B_damagedvar_output
    []
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 40
        show = 'Cd_aux'
    [../]
[]