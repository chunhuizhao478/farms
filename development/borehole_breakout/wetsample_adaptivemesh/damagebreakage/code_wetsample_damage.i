[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../../meshfile/mesh_adaptive.msh'
    []
    displacements = 'disp_x disp_y disp_z'
[]
    
[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
    PorousFlowDictator = dictator #All porous modules must contain

    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 15.62e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 19.92e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.9
    
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
    Cd_constant = 60

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
    C_g = 1e-12 #
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    #coefficient of energy ratio Fb/Fs = chi < 1
    chi = 0.8
    
    #
    D = 0    
[]
  
[Variables]
    #displacement components
    [disp_x]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    [disp_z]
        order = FIRST
        family = LAGRANGE
        scaling = 1E-6
    []
    #pore pressure
    [pp]
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    [] 
    #
    [./biot_modulus_aux]
        order = CONSTANT
        family = MONOMIAL      
    []   
    [./effective_perm00_aux]
        order = CONSTANT
        family = MONOMIAL      
    [] 
    [./effective_perm11_aux]
        order = CONSTANT
        family = MONOMIAL      
    []
    [./effective_perm22_aux]
        order = CONSTANT
        family = MONOMIAL      
    []       
[]
    
#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-2.6477e-5 - 3.3e-7 * t'
    []
    #strain
    [func_strain_xx]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_00
    [../]
    [func_strain_xy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_01
    [../]
    [func_strain_xz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_02
    [../]
    [func_strain_yy]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_11
    [../]
    [func_strain_yz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_12
    [../]
    [func_strain_zz]
        type = SolutionFunction
        solution = init_sol_components
        from_variable = elastic_strain_22
    [../]      
[]
  
[Kernels]
    #effective stress tensor
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        use_displaced_mesh = false
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        use_displaced_mesh = false
    []
    #pressure coupling on stress tensor
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_y
        component = 1
    []
    [poro_z]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_z
        component = 2
    []
    #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        multiply_by_density = false
        biot_coefficient = 0.5
        coupling_type = HydroMechanical
        variable = pp
    []
    #flux * grad(test)
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        multiply_by_density = false
        variable = pp
        gravity = '0 0 0'
    []
[]
  
[AuxKernels]
    [vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_END'
    []
    #
    [biot_modulus]
        type = MaterialRealAux
        property = PorousFlow_constant_biot_modulus_qp
        variable = biot_modulus_aux
        execute_on = 'TIMESTEP_END'
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
    #applied load on inner boundary pore pressure
    [applied_pore_pressure]
        type = DirichletBC
        variable = pp
        boundary = 5
        value = 3.4e6
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
          displacements = 'disp_x disp_y disp_z'
        [../]
        [./inner_boundary]
          boundary = 5
          factor = 3.4e6
          displacements = 'disp_x disp_y disp_z'
        [../]
    []
[]
    
[Materials]
    #damage material
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBMDebug
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 xi stress'
        block = '3'
        outputs = exodus
        ####add porous flow coupling####
        porous_flow_coupling = true
        crack_surface_roughness_correction_factor = 1.0
        length_scale = 0.001
        intrinsic_permeability = 1E-17
    []
    [dummy_matprop]
        type = GenericConstantMaterial
        prop_names = 'initial_damage initial_breakage shear_stress_perturbation'
        prop_values = '0.0 0.0 0.0'  
    [] 
    #elastic material
    [stress_elastic]
        type = ComputeLinearElasticStress
        block = '1 2'
        output_properties = 'elastic_strain stress'
        outputs = exodus
    [] 
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        eigenstrain_names = static_initial_strain_tensor
    []
    [./elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 48.5e9
        poissons_ratio = 0.22
    [../]
    [./initial_strain]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_strain_tensor
        tensor_functions = 'func_strain_xx     func_strain_xy      func_strain_xz 
                            func_strain_xy     func_strain_yy      func_strain_yz
                            func_strain_xz     func_strain_yz      func_strain_zz'
    []
    [density]
        type = GenericConstantMaterial
        prop_names = 'density'
        prop_values = '2640'
    []
    #
    [temperature]
        type = PorousFlowTemperature
    []
    [eff_fluid_pressure_qp]
        type = PorousFlowEffectiveFluidPressure
    []
    #compute volumetric strain and its rate
    [vol_strain]
        type = PorousFlowVolumetricStrain
        outputs = exodus
    []
    #This Material is used for the fully saturated single-phase situation "
    #"where porepressure is the primary variable", saturation = 1.0
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = pp
    []
    #List of variables that represent the mass fractions.
    #If no "variables are provided then num_phases=1=num_components."
    [massfrac]
        type = PorousFlowMassFraction
    []
    #compute porosity
    [porosity]
        type = PorousFlowPorosityConst # only the initial value of this is ever used
        porosity = 0.008 #slide
    []
    #########################################################
    ##### Compute Damaged Permeability and Biot Modulus #####
    #########################################################
    #compute permeability
    [permeability]
        type = FarmsPorousFlowPermeabilityDamaged
        block = '3'
    []
    #compute biot modulus #include damaged solid compliance
    [biot_modulus]
        type = FarmsPorousFlowDamagedBiotModulus
        biot_coefficient = 0.5
        solid_bulk_compliance = 3.46e-11 #calculated
        fluid_bulk_modulus = 2.0e+9
        block = '3'
    []
    ##########################################################
    ##### Compute Constant Permeability and Biot Modulus #####
    ##########################################################
    #comopute permeability
    [permeability_constant]
        type = PorousFlowPermeabilityConst
        permeability = '1E-17 0 0 0 1E-17 0 0 0 1E-17' #slide
        block = '1 2'
    []
    #compute biot modulus
    [biot_modulus_constant]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.5 #paper
        solid_bulk_compliance = 3.46e-11 #calculated
        fluid_bulk_modulus = 2.2e+9
        block = '1 2'
    []    
    ##########################################################
    #Compute density and viscosity
    [simple_fluid_qp]
        type = PorousFlowSingleComponentFluid
        fp = the_simple_fluid
        phase = 0
    []
    #define initial bulk modulus material property (needed by ElkComputeSmearedCrackingStress)
    [solid_bulk_modulus_compliance]
        type = GenericConstantMaterial
        prop_names = solid_bulk_modulus_compliance
        prop_values = 3.46e-11 #calculated
    []
    #define relative permeability as 1 (used in PorousFlowDarcyVelocityComponent)
    [relperm]
        type = PorousFlowRelativePermeabilityConst
        phase = 0
        kr = 1
    []
[]
  
#provide fluid properties for porous flow 
[FluidProperties]
    [the_simple_fluid]
        type = SimpleFluidProperties
        bulk_modulus = 2.2e+9
        density0 = 1000
        thermal_expansion = 0
        viscosity = 1e-3
    []
[]
  
[Preconditioning]
    [smp]
        type = SMP
        full = true
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
    [pp_ic]
        type = SolutionIC
        variable = pp
        solution_uo = init_sol_components
        from_variable = pp
    []
[]
  
#this user object must contain for porous flow
[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pp'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
    [./init_sol_components]
        type = SolutionUserObject
        mesh = ./code_wetsample_static_out.e
        system_variables = 'disp_x disp_y disp_z pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22'
        timestep = LATEST
        force_preaux = true
    [../]
[]
    
[Executioner]
    type = Transient
    solve_type = Newton
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -sub_pc_factor_shift_type'
    # petsc_options_value = 'lu       mumps NONZERO'
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # automatic_scaling = true
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 30
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0
    end_time = 4000
    dt = 0.5
    [./TimeIntegrator]
        type = ImplicitEuler
    [../]
[]

[Outputs] 
    exodus = true
    time_step_interval = 50
    show = 'stress_22 B alpha_damagedvar xi eps_e_22 vel_x vel_y vel_z pp biot_modulus_aux'
    [./csv]
        type = CSV
        time_step_interval = 1
        show = 'strain_z react_z'
    [../]
    [out]
        type = Checkpoint
        time_step_interval = 50
        num_files = 2
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