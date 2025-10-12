[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../../meshfile/mesh_adaptive.msh'
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
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 70

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
    beta_width = 0.05  #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
    
    #
    D = 0


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

    anand_param_go_mat = 0.25
    anand_param_eta_cv_mat = 0.01
    anand_param_p_mat = 1
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
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
    [get_xi]
        type = MaterialRealAux
        variable = xi_aux
        property = xi
        block = '3'
    []
    [get_I2]
        type = MaterialRealAux
        variable = I2_aux
        property = I2
        block = '3'
    [] 
[]

[Kernels]
    [dispkernel_x]
        type = TotalStressDivergenceTensor
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        block = '3'
    []
    [dispkernel_y]
        type = TotalStressDivergenceTensor
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        block = '3'
    []
    [dispkernel_z]
        type = TotalStressDivergenceTensor
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        block = '3'
    []
    [./mass1]
        type = SmallStrainFluidSolidCoupling
        variable = porepressure
        block = '3'
    [../]
    [./mass2]
        type = SmallStrainPorePressureTimeDerivative
        variable = porepressure
        block = '3'
    [../]
    [./darcy_flow]
        type = SmallStrainFluidDiffusion
        variable = porepressure
        block = '3'
    []
    [./darcy_flow_granular]
        type = SmallStrainFluidDiffusionGranular
        variable = porepressure
        block = '3'
    []
    [./plastic_volumetric]
        type = SmallStrainPlasticVolumetricStrainCoupling
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
    [./grad_stress_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        block = '1 2'
  [../]
  [./grad_stress_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        block = '1 2'
  [../]
  [./grad_stress_z]
        type = StressDivergenceTensors
        variable = disp_z
        component = 2
        block = '1 2'
  [../]
  [./poro_x]
        type = PoroMechanicsCoupling
        variable = disp_x
        component = 0
        block = '1 2'    
  [../]
  [./poro_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        component = 1
        block = '1 2'    
  [../]
  [./poro_z]
        type = PoroMechanicsCoupling
        variable = disp_z
        component = 2
        block = '1 2'    
  [../]
  [./poro_timederiv]
        type = PoroFullSatTimeDerivative
        variable = porepressure
        block = '1 2' 
  [../]
  [./darcy_flow2]
        type = CoefDiffusion
        variable = porepressure
        coef = 1E-17
        block = '1 2' 
  [../]
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2640'
    []
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        # outputs = exodus
    [] 
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBM
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 stress'
        block = '3'
        outputs = exodus
    [] 
    [stress_elastic]
        type = ComputeLinearElasticStress
        block = '1 2'
        output_properties = 'elastic_strain stress'
        outputs = exodus
    []
    [elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        block = '1 2'
        youngs_modulus = 48.5e9
        poissons_ratio = 0.22
    []
    [./poro_material]
        type = PoroFullSatMaterial
        block = '1 2'
        porosity0 = 0.008
        biot_coefficient = 0.4264
        solid_bulk_compliance = 1.984914649E-11
        fluid_bulk_compliance = 4.545454545E-10
        constant_porosity = true
    [../]
    [dummy_matprop]
        type = GenericConstantMaterial
        prop_names = 'initial_damage initial_breakage shear_stress_perturbation'
        prop_values = '0.0 0.0 0.0'  
    []    
[]  

#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-3.3e-5 - 3.3e-7 * t'
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
    end_time = 4000 #extend the time
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 30
    nl_abs_tol = 1e-8
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'none'
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1
        cutback_factor_at_failure = 0.5
        optimal_iterations = 10
        growth_factor = 1.25
        max_time_step_bound = 3
    []
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
[]

[Outputs] 
    exodus = true
    time_step_interval = 5
    show = 'stress_22 porepressure B alpha_damagedvar xi eps_e_22 vel_x vel_y vel_z'
    [./csv]
        type = CSV
        time_step_interval = 1
        show = 'strain_z react_z'
    [../]
    [out]
        type = Checkpoint
        time_step_interval = 20
        num_files = 2
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
    [fix_top_x]
        type = DirichletBC
        variable = disp_x
        boundary = 6
        value = 0
    []
    [fix_top_y]
        type = DirichletBC
        variable = disp_y
        boundary = 6
        value = 0
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

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_main_out.e'
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
      variable = xi_aux
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
      stress_tensor = stress
      boundary = 6
    [../]
    [./strain_z]
        type = FunctionValuePostprocessor
        function = applied_load_top
    []
[]
