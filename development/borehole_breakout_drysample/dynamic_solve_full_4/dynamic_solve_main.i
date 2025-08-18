#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_adaptive.msh'
        show_info = true
    [] 

[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    porepressure = 'porepressure'
      
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 15.62e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 19.92e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8073
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.8073
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8

    # Water bulk modulus (2.2 GPa)
    fluid_bulk_modulus = 2.2e9      

     # Initial permeability (1 milli-darcy) 
    permeability_solid_o = 1e-20 

     # Initial porosity (15%)
    porosity_solid_o = 0.008   
     
     # Granular bulk modulus (15 GPa)      
    solid_bulk_modulus_g = 50.38e9 

    # Solid grains bulk modulus (36 GPa - typical for quartz)  
    solid_bulk_modulus_s = 50.38e9 

    permeability_evolution_with_damage = 3
    initial_grain_size = 1.3
    ultimate_grain_size = 0.25
    initial_viscosity_fluid = 1e-3
    

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
    [porepressure]
        order = FIRST
        family = LAGRANGE
    []
    
[]

[AuxVariables]
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
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
    #
    [alpha_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    [B_damagedvar_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
    [I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [deviatroic_strain_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [structural_stress_coefficient_aux]
        order = FIRST
        family = MONOMIAL
    []
    #
    [gradx_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    [grady_alpha_damagedvar]
        order = CONSTANT
        family = MONOMIAL
    []
    #spatial damage parameters
    [cg_aux]
        order = FIRST
        family = LAGRANGE
    []
    #
    [nonlocal_xi]
        order = FIRST
        family = MONOMIAL
    []
    
[]

[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []  
    #
    [get_xi]
        type = MaterialRealAux
        variable = xi_aux
        property = strain_invariant_ratio
        block = '3'
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_aux
        property = second_elastic_strain_invariant
        block = '3'
    [] 
    [get_deviatroic_strain_rate]
        type = MaterialRealAux
        variable = deviatroic_strain_rate_aux
        property = deviatroic_strain_rate
        block = '3'
    []
    #
    [get_nonlocal_xi]
        type = MaterialRealAux
        variable = nonlocal_xi
        property = eqstrain_nonlocal
    []
[]

[Kernels]
    [grad_stress_x]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_x
        component = 0
        large_kinematics = true
    []
    [grad_stress_y]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_y
        component = 1
        large_kinematics = true
    []
    [grad_stress_z]
        type = TotalLagrangianTotalStressDivergence
        variable = disp_z
        component = 2
        large_kinematics = true
    []
    [./mass1]
        type = FluidSolidCoupling
        variable = porepressure
    [../]
    [./mass2]
        type = PorePressureTimeDerivative
        variable = porepressure
    [../]
    [./darcy_flow]
        type = FluidDiffusion
        variable = porepressure
        large_kinematics = false
    []
    [./darcy_flow_granular]
        type = FluidDiffusionGranular
        variable = porepressure
        block = '3'
    []
    [./plastic_volumetric]
        type = PlasticVolumetricStrainCoupling
        variable = porepressure
        block = '3'
    []
    [./inertia_x]
        type = InertialForce
        variable = disp_x
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
    [./inertia_y]
        type = InertialForce
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
        eta = 0
    []    
    [./inertia_z]
        type = InertialForce
        variable = disp_z
        acceleration = accel_z
        velocity = vel_z
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2640'
    []
    [strain]
        type = ComputeLagrangianStrain
        large_kinematics = true
        # outputs = exodus
    []
    #
    [damage_mat]
        type = DiffusedDamageBreakageMaterialMainApp
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        structural_stress_coefficient = structural_stress_coefficient_aux
        #build L matrix using velocity
        vel_x = vel_x
        vel_y = vel_y
        vel_z = vel_z
        block = '3'
    [] 
    [stress_medium]
        type = ComputeLagrangianDamageBreakageStressPK2Diffused
        large_kinematics = true
        output_properties = 'pk2_stress green_lagrange_elastic_strain plastic_strain total_lagrange_strain strain_invariant_ratio'
        outputs = exodus
        block = '3'
    []
    #elastic material
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 48.5e9
        poissons_ratio = 0.22
    []
    [porous_prop]
        type =   IntactPorousSolidProperties
        block = '1 2'
    []
    [compute_stress]
        type = ComputePoroStVenantKirchhoffStress
        large_kinematics = true
        output_properties = 'green_lagrange_strain pk2_stress'
        outputs = exodus
        block = '1 2'
    []
    #strain invariant ratio
    [comp_strain_invariant_ratio]
        type = ComputeXi 
        output_properties = 'strain_invariant_ratio'
        outputs = exodus
        block = '1 2'
    []
    #nonlocal eqstrain
    [nonlocal_eqstrain]
        type = ElkNonlocalEqstrain
        average_UO = eqstrain_averaging
        output_properties = 'eqstrain_nonlocal'
        outputs = exodus
    []
    #shear stress perturbation
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'shear_stress_perturbation damage_perturbation'
        prop_values = '0.0 0.0'
    []
[] 

[UserObjects]
    [eqstrain_averaging] #length scale = radius = grain size 
        type = ElkRadialAverage
        length_scale = 0.0013
        prop_name = strain_invariant_ratio
        radius = 0.0013
        weights = BAZANT
        execute_on = LINEAR
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
    start_time = -1e-12
    end_time = 1e10
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 30
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  boomeramg True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    nl_forced_its = 3
    line_search = 'bt'
    # dt = 10
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1
        cutback_factor_at_failure = 0.75
        optimal_iterations = 25
        growth_factor = 1.25
        max_time_step_bound = 10
    []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 3 ###
        show = 'disp_z porepressure vel_x vel_y vel_z alpha_damagedvar_aux B_damagedvar_aux xi_aux deviatroic_strain_rate_aux nonlocal_xi pk2_stress_22 green_lagrange_elastic_strain_22 plastic_strain_22 total_lagrange_strain_22'
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
    [../]
    [checkpoint]
        type = Checkpoint
        time_step_interval = 20
        num_files = 2
    []
[]

#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-3.3e-5 - 3.3e-7 * t'
    []
[]


[BCs]
    #fix bottom boundary
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 7
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 7
        value = 0
    []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = 7
        value = 0
    []
    #applied load on top boundary
    [applied_top_z_dispload]
        type = FunctionDirichletBC
        variable = disp_z
        boundary = 6
        function = applied_load_top
    [] 
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = 4
          factor = 20.6e6
          displacements = 'disp_x disp_y'
        [../]
    []
    [./PorePressure]
          type = FunctionDirichletBC
          boundary = 5
          variable = porepressure
          function = 3.4e6
    []
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'dynamic_solve_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
        sub_cycling = true
        clone_parent_mesh = true
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_damagedvar_sub B_damagedvar_sub structural_stress_coefficient_sub'
        variable = 'alpha_damagedvar_aux B_damagedvar_aux structural_stress_coefficient_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'I2_aux nonlocal_xi deviatroic_strain_rate_aux'
        variable = 'I2_sub_aux xi_sub_aux deviatroic_strain_rate_sub_aux'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_exodus.e'
      system_variables = 'disp_x disp_y disp_z porepressure initial_xi_aux initial_I2_aux'
      timestep = LATEST
      force_preaux = true
      execute_on = 'INITIAL'
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
    [strain_invariant_ratio_ic]
      type = SolutionIC
      variable = nonlocal_xi
      solution_uo = init_sol_components
      from_variable = initial_xi_aux
    []
    [I2_aux_ic]
      type = SolutionIC
      variable = I2_aux
      solution_uo = init_sol_components
      from_variable = initial_I2_aux
    []  
    [pore_pressure_aux_ic]
      type = SolutionIC
      variable = porepressure
      solution_uo = init_sol_components
      from_variable = porepressure
    []  
[]

#compute the reaction force on the top boundary
[Postprocessors]
    [./react_z]
      type = SidesetReaction
      direction = '0 0 1'
      stress_tensor = pk2_stress
      boundary = 6
    [../]
    [./strain_z]
        type = FunctionValuePostprocessor
        function = applied_load_top
    []
[]