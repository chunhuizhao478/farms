#solid properties
#-------------------------------------------------#
solid_density = 2640
youngs_modulus = 48.5e9
poissons_ratio = 0.22
solid_bulk_compliance = 3.46e-11
lambda_o = ${fparse youngs_modulus*poissons_ratio/(1+poissons_ratio)/(1-2*poissons_ratio)}
shear_modulus_o = ${fparse youngs_modulus/(2*(1+poissons_ratio))}
length_scale = 3e-3
crack_surface_roughness_correction_factor = 1.0
#-------------------------------------------------#

#damage-breakage properties
#-------------------------------------------------#
xi_o = -0.8073
xi_d = -0.8073
Cg = 1e-12
#-------------------------------------------------#

#fluid properties
#-------------------------------------------------#
porosity = 0.008
permeability = '1E-20 0 0 0 1E-20 0 0 0 1E-20'
intrinsic_permeability = 1E-20
fluid_density = 1000
viscosity = 1e-3
biot_coefficient = 0.5
fluid_bulk_modulus = 2.2e+9
#-------------------------------------------------#

#boundary conditions
#-------------------------------------------------#
outer_confinement_pressure = 20.6e6
inner_confinement_pressure = 3.4e6
#-------------------------------------------------#

#implicit continuum damage-breakage model dynamics
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_adaptive_test.msh'
    [] 
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    PorousFlowDictator = dictator #All porous modules must contain
      
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = ${lambda_o}
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = ${shear_modulus_o}
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = ${xi_o}
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = ${xi_d}
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = ${Cg}
    
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
    [dispkernel_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        variable = disp_z
        component = 2
    []
    #pressure coupling on stress tensor
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = ${biot_coefficient}
        variable = disp_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = ${biot_coefficient}
        variable = disp_y
        component = 1
    []
    [poro_z]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = ${biot_coefficient}
        variable = disp_z
        component = 2
    []
    #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        multiply_by_density = false
        biot_coefficient = ${biot_coefficient}
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

[Materials]
    [strain]
        type = ComputeSmallStrain
        eigenstrain_names = static_initial_strain_tensor
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBMDiffused
        alpha_damagedvar_aux = alpha_damagedvar_aux
        B_damagedvar_aux = B_damagedvar_aux
        output_properties = 'stress elastic_strain_tensor plastic_strain_tensor total_strain_tensor strain_invariant_ratio'
        outputs = exodus
        block = '3'
        #porous flow coupling
        porous_flow_coupling = true
        crack_surface_roughness_correction_factor = ${crack_surface_roughness_correction_factor}
        length_scale = ${length_scale}
        intrinsic_permeability = ${intrinsic_permeability}
    []
    #elastic material
    [elastic_tensor]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = ${youngs_modulus}
        poissons_ratio = ${poissons_ratio}
    []
    [compute_stress]
        type = ComputeLinearElasticStress
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
    [damage_perturbation]
        type = PerturbationRadial
        nucl_center = '0 0 0'
        peak_value = 0
        thickness = 200
        length = 2000
        duration = 1.0
        perturbation_type = 'shear_stress'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
        outputs = exodus
    []
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
        prop_values = ${solid_density}
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
        porosity = ${porosity}
    [] 
    #########################################################
    ##### Compute Damaged Permeability and Biot Modulus #####
    #########################################################
    # #compute permeability
    # [permeability]
    #     type = FarmsPorousFlowPermeabilityDamaged
    #     block = '3'
    # []
    # #compute biot modulus #include damaged solid compliance
    # [biot_modulus]
    #     type = FarmsPorousFlowDamagedBiotModulus
    #     biot_coefficient = 0.5
    #     solid_bulk_compliance = 3.46e-11 #calculated
    #     fluid_bulk_modulus = 2.2e+9
    #     block = '3'
    # []
    ##########################################################
    ##### Compute Constant Permeability and Biot Modulus #####
    ##########################################################
    #compute permeability
    [permeability_constant]
        type = PorousFlowPermeabilityConst
        permeability = ${permeability}
        block = '1 2 3'
    []
    #compute biot modulus
    [biot_modulus_constant]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = ${biot_coefficient}
        solid_bulk_compliance = ${solid_bulk_compliance}
        fluid_bulk_modulus = ${fluid_bulk_modulus}
        block = '1 2 3'
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
        prop_values = ${solid_bulk_compliance}
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
        bulk_modulus = ${fluid_bulk_modulus}
        density0 = ${fluid_density}
        thermal_expansion = 0
        viscosity = ${viscosity}
    []
[]

[UserObjects]
    [eqstrain_averaging] #length scale = radius = grain size 
        type = ElkRadialAverage
        length_scale = ${length_scale}
        prop_name = strain_invariant_ratio
        radius = ${length_scale}
        weights = BAZANT
        execute_on = TIMESTEP_END
    []
[]

#18.2e6 * 0.1 / 48.5e9 = 3.7525e-5 applied displacement (seating load)
[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = 'if (t > 1e-3, -2.6477e-5 - 3.3e-7 * t, -2.6477e-5)'
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

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Controls] # turns off inertial terms for the SECOND time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/mass0'
    start_time = -1e-12
    end_time = 1e-3 # dt used in the simulation
  []
[../]
  
[Executioner]
    type = Transient
    # solve_type = 'NEWTON'
    solve_type = 'PJFNK'
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
    # automatic_scaling = true
    # nl_forced_its = 3
    line_search = 'bt'
    # dt = 10
    verbose = true
    [TimeStepper]
        type = FarmsIterationAdaptiveDT
        dt = 1e-3
        cutback_factor_at_failure = 0.5
        optimal_iterations = 20
        growth_factor = 1.25
        max_time_step_bound = 100
    []
    [./TimeIntegrator]
        type = ImplicitEuler
        # type = BDF2
        # type = CrankNicolson
    [../]
[]

[Outputs]
    [./exodus]
        type = Exodus
        time_step_interval = 1 ###
        # show = 'vel_x vel_y vel_z alpha_damagedvar_aux B_damagedvar_aux xi_aux deviatroic_strain_rate_aux nonlocal_xi stress_22 elastic_strain_tensor_22 plastic_strain_tensor_22 total_strain_tensor_22'
    [../]
    [./csv]
        type = CSV
        time_step_interval = 1
    [../]
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
          factor = ${outer_confinement_pressure}
          displacements = 'disp_x disp_y'
        [../]
        [./inner_boundary]
          boundary = 5
          factor = ${inner_confinement_pressure}
          displacements = 'disp_x disp_y disp_z'
        [../]
    []
    #inner pressure
    # [inner_pressure]
    #     type = PorousFlowPiecewiseLinearSink
    #     variable = pp
    #     boundary = '5'
    #     pt_vals = '-1e9 1e9' # x coordinates defining g
    #     multipliers = '-1e9 1e9' # y coordinates defining g
    #     PT_shift = ${inner_confinement_pressure}
    #     flux_function = 1E-6 # Variable C
    #     fluid_phase = 0
    # []
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
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pp'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
    [./init_sol_components]
        type = SolutionUserObject
        mesh = ../static_solve/static_solve_out.e
        system_variables = 'disp_x disp_y disp_z pp elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22 initial_I2_aux initial_xi_aux'
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
    [pp_ic]
        type = SolutionIC
        variable = pp
        solution_uo = init_sol_components
        from_variable = pp
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

[VectorPostprocessors]
    [output_vel_x]
        type = NodalValueSampler
        variable = vel_x
        sort_by = 'id'
    []
    [output_vel_y]
        type = NodalValueSampler
        variable = vel_y
        sort_by = 'id'
    []
    [output_vel_z]
        type = NodalValueSampler
        variable = vel_z
        sort_by = 'id'
    []
[]