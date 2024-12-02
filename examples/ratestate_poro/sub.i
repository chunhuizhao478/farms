[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = './Inclined_fault_with_injection.msh'
    []
    [subdomain1]
        input = msh
        type = SubdomainBoundingBoxGenerator
        bottom_left = '-5000 -5000 0'
        top_right = '5000 5000 0'
        block_id = 0
    []
    [./inner_block]
        type = ParsedSubdomainMeshGenerator
        input = subdomain1
        combinatorial_geometry = 'x > -1732.05 & x < 1732.05'
        block_id = 1
    []
    [./fault_block_upper]
        type = ParsedSubdomainMeshGenerator
        input = inner_block
        combinatorial_geometry = 'y > (-0.577350269 * x) & x > -1732.05 & x < 1732.05'
        block_id = 2
    []
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = fault_block_upper
        split_interface = true
        add_interface_on_two_sides = true
        block_pairs = '1 2'
    []
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y' 

    PorousFlowDictator = dictator

    q = 0.0
    elem_size = 4

    
    ##element length (m)
    len = 4
    
    ##rate-and-state coefficients
    f_o = 0.6
    rsf_a = 0.005
    rsf_b = 0.02
    rsf_L = 0.0002
    delta_o = 1e-9

    ##initial normal traction (Pa)
    Tn_o = 0 # check if we should include the far field tectonic loading

    ##initial shear traction (Pa)
    Ts_o = 0 # check if we should include the far field tectonic loading

    ##initial sliprate (m/s)
    Vini = 2.17e-2 # assumed #  2.17e-2 steady state

    ##initial state variable
    statevarini = 0.00922 #18003426.26 #  steady state
    # statevarini = 4068457.682 # if sinitial sliprate is 1e-9
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'p_sub'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = domain_Block2
      execute_on = 'initial TIMESTEP_BEGIN'
   [../]
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0.0
    bulk_modulus = 2.5e9
    viscosity = 0.001
    density0 = 1000
  []
[]

[Variables]
    [./disp_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p_sub]
        order = FIRST
        family = LAGRANGE
    [../]
[]


[AuxVariables]
    [./nodal_area]
        order = SECOND
        family = LAGRANGE
    [../]
    #fault displacement residual from interfacekernel
    [disp_plusminus_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [disp_plusminus_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    #fault displacement residual from interfacekernel
    [fluid_disp_plusminus_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_disp_plusminus_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_vel_plusminus_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_vel_plusminus_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plusminus_sub_scaled_x]
        order = SECOND
        family = LAGRANGE
    []
    [disp_plusminus_sub_scaled_y]
        order = SECOND
        family = LAGRANGE
    []
    #side element volume 
    [element_side_volume]
        order = SECOND
        family = LAGRANGE
    []
    #residual received from mainApp
    [resid_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [resid_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    [./resid_visco_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_visco_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    #restoration pressure force (tag after solve)
    [./resid_pressure_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    #residual received from mainApp (damping)
    [resid_damp_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [resid_damp_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    [./resid_pressdamp_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressdamp_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = SECOND
        family = MONOMIAL
    []
    #traction
    [traction_sub_strike]
        order = SECOND
        family = MONOMIAL
    []
    [traction_sub_normal]
        order = SECOND
        family = MONOMIAL
    []
    #sliprate
    [sliprate_sub_strike]
        order = SECOND
        family = MONOMIAL
    []
    [sliprate_sub_normal]
        order = SECOND
        family = MONOMIAL
    []
    #slip
    [slip_sub_strike]
        order = SECOND
        family = MONOMIAL
    []
    [slip_sub_normal]
        order = SECOND
        family = MONOMIAL
    []
    #statevar
    [statevar_sub]
        order = SECOND
        family = MONOMIAL
    []
    [p_sub_main]
        order = FIRST
        family = LAGRANGE
    []
[]



[AuxKernels]
    #retrieve displacement vector by scaling the residuals
    [scale_disp_residual_x]
        type = ScaleVarAux
        variable = disp_plusminus_sub_scaled_x
        coupled = disp_plusminus_sub_x
        scale = element_side_volume
        execute_on = 'TIMESTEP_END'
    []
    [scale_disp_residual_y]
        type = ScaleVarAux
        variable = disp_plusminus_sub_scaled_y
        coupled = disp_plusminus_sub_y
        scale = element_side_volume
        execute_on = 'TIMESTEP_END'
    []
    #const element_side_volume
    [const_element_side_volume]
        type = ConstantAux
        variable = element_side_volume
        value = 100
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #retrieve fault displacement residual vector using tagging
    [restore_x]
        type = TagVectorAux
        vector_tag = 'restore_tag_x'
        v = 'disp_sub_x'
        variable = 'disp_plusminus_sub_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_y]
        type = TagVectorAux
        vector_tag = 'restore_tag_y'
        v = 'disp_sub_y'
        variable = 'disp_plusminus_sub_y'
        execute_on = 'TIMESTEP_END'
    []
    ##scale the displacement by surface area

    ##initial shear traction
    [StrikeShearStress]
        type = FunctionAux
        variable = ini_shear_stress_perturb
        function = func_initial_strike_shear_stress
        execute_on = 'TIMESTEP_BEGIN'
    []
    ##output
    [output_traction_strike]
        type = MaterialRealAux
        property = traction_strike
        variable = traction_sub_strike
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_traction_normal]
        type = MaterialRealAux
        property = traction_normal
        variable = traction_sub_normal
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_sliprate_strike]
        type = MaterialRealAux
        property = sliprate_strike
        variable = sliprate_sub_strike
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_sliprate_normal]
        type = MaterialRealAux
        property = sliprate_normal
        variable = sliprate_sub_normal
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_strike]
        type = MaterialRealAux
        property = slip_strike
        variable = slip_sub_strike
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_normal]
        type = MaterialRealAux
        property = slip_normal
        variable = slip_sub_normal
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
    [output_statevar]
        type = MaterialRealAux
        property = statevar
        variable = statevar_sub
        check_boundary_restricted = false
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
    []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'domain_Block2'
        strain = SMALL
    [../]
[]

[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_sub_x
        component = 0
        displacements = 'disp_sub_x disp_sub_y'
        use_displaced_mesh = false    
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_sub_y
        component = 1
        displacements = 'disp_sub_x disp_sub_y'
        use_displaced_mesh = false    
    [../]
    [viscoelastic_x]
        type = ViscoelasticStressKernel
        variable = disp_sub_x
        component = 0
        kelvin_voigt_viscosity = 10  
        shear_modulus = 1e9        
    []
    [viscoelastic_y]
        type = ViscoelasticStressKernel
        variable = disp_sub_y
        component = 1
        kelvin_voigt_viscosity = 10  
        shear_modulus = 1e9        
    []  
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_sub_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_sub_y
        use_displaced_mesh = false
    [../]
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_sub_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_sub_y
        component = 1
    []
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        biot_coefficient = 0.5
        coupling_type = HydroMechanical
        variable = p_sub
        multiply_by_density = false
    []
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p_sub
        gravity = '0 0 0'
        multiply_by_density = false
    []
    [./Reactionx]
        type = StiffPropDamping
        variable = disp_sub_x
        component = 0
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = disp_sub_y
        component = 1
    []
    [./Reactionx_P]
        type = PoroPropDampingCoupling
        variable = disp_sub_x
        component = 0
        porepressure = p_sub
    []
    [./Reactiony_p]
        type = PoroPropDampingCoupling
        variable = disp_sub_y
        component = 1
        porepressure = p_sub
    []
    
[]


[Problem]
    extra_tag_vectors = 'restore_tag_x restore_tag_y'
[]


[InterfaceKernels]
    #apply displacement prediction and retrieve its residuals
    [./ratestate_x]
        type = RateStateInterfaceKernelGlobalx
        variable = disp_sub_x
        neighbor_var = disp_sub_x
        extra_vector_tags = 'restore_tag_x'
        boundary = 'domain_Block2'
        y_var = disp_sub_y
    []
    [./ratestate_y]
        type = RateStateInterfaceKernelGlobaly
        variable = disp_sub_y
        neighbor_var = disp_sub_y
        extra_vector_tags = 'restore_tag_y'
        boundary = 'domain_Block2'
        x_var = disp_sub_x
    []
[]


[Materials]
    [temperature]
        type = PorousFlowTemperature
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 20e9
        poissons_ratio = 0.25
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [Strain]
        type = ComputeSmallStrain
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2350
    []
    [rho_f]
        type = GenericConstantMaterial
        prop_names = rhof
        prop_values = 1000
    []
    [biot_coeff]
        type = GenericConstantMaterial
        prop_names = 'biot_coefficient'
        prop_values = 0.5    
        block = '0 1 2'        
    []
    [eff_fluid_pressure_qp]
        type = PorousFlowEffectiveFluidPressure
    []
    [vol_strain]
        type = PorousFlowVolumetricStrain
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = p_sub
    []
    [massfrac]
        type = PorousFlowMassFraction
    []
    [simple_fluid_qp]
        type = PorousFlowSingleComponentFluid
        fp = simple_fluid
        phase = 0
    []
    [porosity]
        type = PorousFlowPorosityConst # only the initial value of this is ever used
        porosity = 0.1
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.5
        solid_bulk_compliance = 3.75e-11  
        fluid_bulk_modulus = 2.5e9        
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '1.128052989e-12 0 0   0 1.128052989e-12 0   0 0 1.128052989e-12'
    []
    [./czm_mat]
        type = PoroRateStateFrictionLaw2DAsBC
        reaction_rsf_x  = resid_sub_x
        reaction_rsf_y  = resid_sub_y
        reaction_visco_rsf_x = resid_visco_sub_x
        reaction_visco_rsf_y = resid_visco_sub_y
        reaction_rsf_pressure_x = resid_pressure_sub_x
        reaction_rsf_pressure_y = resid_pressure_sub_y
        reaction_damp_x = resid_damp_sub_x
        reaction_damp_y = resid_damp_sub_y
        reaction_pressdamp_x = resid_pressdamp_sub_x
        reaction_pressdamp_y = resid_pressdamp_sub_y
        fluid_disp_x = fluid_disp_plusminus_sub_x
        fluid_disp_y = fluid_disp_plusminus_sub_y
        fluid_vel_x = fluid_vel_plusminus_sub_x
        fluid_vel_y = fluid_vel_plusminus_sub_y
        nodal_area = nodal_area
        interface_pressure = p_sub_main
        permeability_type = 'Semi_permeable'
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'domain_Block2'
        output_properties = 'sliprate_strike slip_strike statevar traction_strike traction_normal alongfaultdisp_strike_plus alongfaultdisp_strike_minus'
        # outputs = exodus
    [../]
[]

[UserObjects]
    #compute element side volume (using CONTACT modulus)
    # [element_side_volume]
    #     type = NodalArea
    #     variable = element_side_volume
    #     boundary = 'domain_Block2 domain_Block2'
    #     execute_on = 'initial TIMESTEP_BEGIN'
    # []
    [recompute_residual_tag_x]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag_x'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_y]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag_y'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
[]

[Functions]
    [./func_initial_strike_shear_stress]
        type = InitialStrikeShearStressPerturbRSF2D
    []
[]

#[Preconditioning]
#  [./smp]
#    type = SMP
#    full = true
#  [../]
#[]

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]