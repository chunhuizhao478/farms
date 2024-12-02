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
    displacements = 'disp_x disp_y' 
    PorousFlowDictator = dictator
    q = 0.5
    gravity = '0 0 0'
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
    [./disp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    [./vel_x]
        order = SECOND
        family = LAGRANGE
    []
    [./accel_x]
        order = SECOND
        family = LAGRANGE
    []
    [./vel_y]
        order = SECOND
        family = LAGRANGE
    []
    [./accel_y]
        order = SECOND
        family = LAGRANGE
    []
    [./nodal_area]
        order = SECOND
        family = LAGRANGE
    [../]
    #restoration force (tag after solve)
    [./resid_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_visco_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_visco_y]
        order = SECOND
        family = LAGRANGE
    [../]
    #restoration pressure force (tag after solve)
    [./resid_pressure_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_y]
        order = SECOND
        family = LAGRANGE
    [../] 
    #restoration force for damping (tag after solve)
    [./resid_damp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressdamp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressdamp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plusminus_scaled_x]
        order = SECOND
        family = LAGRANGE
    []
    [disp_plusminus_scaled_y]
        order = SECOND
        family = LAGRANGE
    []
    #traction
    [traction_strike]
        order = SECOND
        family = MONOMIAL
    []
    [traction_normal]
        order = SECOND
        family = MONOMIAL
    []
    #sliprate
    [sliprate_strike]
        order = SECOND
        family = MONOMIAL
    []
    [sliprate_normal]
        order = SECOND
        family = MONOMIAL
    []
    #slip
    [slip_strike]
        order = SECOND
        family = MONOMIAL
    []
    [slip_normal]
        order = SECOND
        family = MONOMIAL
    []
    #statevar
    [statevar]
        order = SECOND
        family = MONOMIAL
    []
    [fluid_disp_plusminus_x]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_disp_plusminus_y]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_vel_plusminus_x]
        order = SECOND
        family = LAGRANGE
    []
    [fluid_vel_plusminus_y]
        order = SECOND
        family = LAGRANGE
    []  
    [./elemental_across_flux_main]
        order = SECOND
        family = MONOMIAL
    []
    [./elemental_across_flux_sec]
        order = SECOND
        family = MONOMIAL
    []
    [./across_flux_main]
        order = SECOND
        family = LAGRANGE
    []
    [./across_flux_sec]
        order = SECOND
        family = LAGRANGE
    []
[]

[AuxKernels]
    #obtain system residuals by tagging
    [restore_x]
        type = TagVectorAux
        vector_tag = 'restore_tag_x'
        v = 'disp_x'
        variable = 'resid_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_y]
        type = TagVectorAux
        vector_tag = 'restore_tag_y'
        v = 'disp_y'
        variable = 'resid_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_visco_x]
        type = TagVectorAux
        vector_tag = 'restore_visco_tag_x'
        v = 'disp_x'
        variable = 'resid_visco_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_visco_y]
        type = TagVectorAux
        vector_tag = 'restore_visco_tag_y'
        v = 'disp_y'
        variable = 'resid_visco_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressurex]
        type = TagVectorAux
        vector_tag = 'restore_pressurex_tag'
        v = 'disp_x'
        variable = 'resid_pressure_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressurey]
        type = TagVectorAux
        vector_tag = 'restore_pressurey_tag'
        v = 'disp_y'
        variable = 'resid_pressure_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampx]
        type = TagVectorAux
        vector_tag = 'restore_dampx_tag'
        v = 'disp_x'
        variable = 'resid_damp_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampy]
        type = TagVectorAux
        vector_tag = 'restore_dampy_tag'
        v = 'disp_y'
        variable = 'resid_damp_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressdampx]
        type = TagVectorAux
        vector_tag = 'restore_pressdampx_tag'
        v = 'disp_x'
        variable = 'resid_pressdamp_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_pressdampy]
        type = TagVectorAux
        vector_tag = 'restore_pressdampy_tag'
        v = 'disp_y'
        variable = 'resid_pressdamp_y'
        execute_on = 'TIMESTEP_END'
    []
    #calc velocity
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
    # Calculate accelerations by taking time derivative of velocities
    [acc_x]
        type = CompAcceleration
        variable = accel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [acc_y]
        type = CompAcceleration
        variable = accel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [flux_main_side]
        type = MaterialRealAux
        property = across_fault_flux
        variable = elemental_across_flux_main
        boundary = 'domain_Block2'
        execute_on = 'TIMESTEP_END'
        check_boundary_restricted = false
    []     
    [flux_sec_side]
        type = MaterialRealAux
        property = across_fault_flux
        variable = elemental_across_flux_sec
        boundary = 'Block2_domain'
        execute_on = 'TIMESTEP_END'
        check_boundary_restricted = false
    [] 
    [nodal_flux_main]
        type = ProjectionAux
        variable = across_flux_main
        v = elemental_across_flux_main
        execute_on = 'TIMESTEP_END'
    []
    [nodal_flux_sec]
        type = ProjectionAux
        variable = across_flux_sec
        v = elemental_across_flux_sec
        execute_on = 'TIMESTEP_END'
    []
[]

[Problem]
    extra_tag_vectors = 'restore_tag_x restore_tag_y restore_visco_tag_x restore_visco_tag_y restore_pressurex_tag restore_pressurey_tag restore_dampx_tag restore_dampy_tag restore_pressdampx_tag restore_pressdampy_tag'
[]

[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false    
        extra_vector_tags = 'restore_tag_y' 
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false    
        extra_vector_tags = 'restore_tag_y' 
    [../]
    [viscoelastic_x]
        type = ViscoelasticStressKernel
        variable = disp_x
        component = 0
        kelvin_voigt_viscosity = 10  
        shear_modulus = 8e9   
        extra_vector_tags = 'restore_visco_tag_x'       
    []
    [viscoelastic_y]
        type = ViscoelasticStressKernel
        variable = disp_y
        component = 1
        kelvin_voigt_viscosity = 10  
        shear_modulus = 8e9     
        extra_vector_tags = 'restore_visco_tag_y'     
    []  
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_y
        use_displaced_mesh = false
    [../]
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.95
        variable = disp_x
        component = 0
        extra_vector_tags = 'restore_pressurex_tag' 
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.95
        variable = disp_y
        component = 1
        extra_vector_tags = 'restore_pressurey_tag' 
    []
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        biot_coefficient = 0.95
        coupling_type = HydroMechanical
        variable = p
        multiply_by_density = false
    []
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p
        multiply_by_density = false
    []
    [./Reactionx]
        type = StiffPropDamping
        variable = disp_x
        component = 0
        extra_vector_tags = 'restore_dampx_tag'
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = disp_y
        component = 1
        extra_vector_tags = 'restore_dampy_tag'
    []
    [./Reactionx_P]
        type = PoroPropDampingCoupling
        variable = disp_x
        component = 0
        porepressure = p
        extra_vector_tags = 'restore_pressdampx_tag'
    []
    [./Reactiony_p]
        type = PoroPropDampingCoupling
        variable = disp_y
        component = 1
        porepressure = p
        extra_vector_tags = 'restore_pressdampy_tag'
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
        prop_values = 0.95    
        block = '0 1 2'        
    []
    [./Transmissibility]
        type = GenericConstantMaterial
        prop_names = Transmissibility
        prop_values = 1e-14
    [../]
    [eff_fluid_pressure_qp]
        type = PorousFlowEffectiveFluidPressure
    []
    [vol_strain]
        type = PorousFlowVolumetricStrain
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = p
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
        biot_coefficient = 0.95
        solid_bulk_compliance = 3.75e-11  # 1/Ks
        fluid_bulk_modulus = 2.5e9        # Kf
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '10e-14 0 0   0 10e-14 0   0 0 10e-14'
    []
    [./flux_main_mat]
        type = SemiPermeableFault
        boundary = 'domain_Block2'
        pressure_secondary = p
        pressure_main = p
    [../]
    [./flux_sec_mat]
        type = SemiPermeableFault
        boundary = 'Block2_domain'
        pressure_secondary = p
        pressure_main =  p
    [../]
   # [fluid_vel_x_aux]
   #     type = PorousFlowDarcyVelocityComponent
   #     variable = fluid_vel_x
   #     component = x
   #     fluid_phase = 0
   # []
[]

[UserObjects]
    [injected_mass]
        type = PorousFlowSumQuantity
    []
    #evalute system residual after system solve before auxkernels
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
    [recompute_residual_visco_tag_x]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_visco_tag_x'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_visco_tag_y]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_visco_tag_y'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
        [recompute_residual_pressure_tag_x]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressurex_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_pressure_tag_y]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressurey_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampx]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampy]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampy_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_pressdampx]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressdampx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_pressdampy]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_pressdampy_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'p'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
    [./nodal_area]
        type = NodalArea
        variable = nodal_area
        boundary = domain_Block2
        execute_on = 'initial TIMESTEP_BEGIN'
    [../]
    [./init_sol_components]
        type = SolutionUserObject
        mesh = 'main_mesh_out.e'  
        system_variables = 'disp_x disp_y p'
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
    [p_ic]
      type = SolutionIC
      variable = p
      solution_uo = init_sol_components
      from_variable = p
    []
[]

[BCs]
    #assign displacement boundary condition
    [./matchval_primary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plusminus_scaled_x
        boundary = 'domain_Block2'
    []
    [./matchval_secondary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plusminus_scaled_x
        boundary = 'Block2_domain'
    []
    [./matchval_primary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_plusminus_scaled_y
        boundary = 'domain_Block2'
    []
    [./matchval_secondary_y]
        type = MatchedValueBC
        variable = disp_y
        v = disp_plusminus_scaled_y
        boundary = 'Block2_domain'
    []
    # Displacement and tractions boundary conditions
    [./disp_bottom]
        type = DirichletBC
        variable = disp_y
        boundary = '1'
        value = 0
    [../]
    [./disp_left]
        type = DirichletBC
        variable = disp_x
        boundary = 'left'
        value = 0
    [../]
    [./traction_right]
        type = NeumannBC
        variable = disp_x
        boundary = 'right'
        value = 0
    [../]
    [./traction_top]
        type = NeumannBC
        variable = disp_y
        boundary = 'top'
        value = 0
    [../]
    [./flux_main]
        type = CoupledVarNeumannBC
        variable = p
        boundary = 'domain_Block2'
        v = across_flux_main
    [../]
    [./flux_sec]
        type = CoupledVarNeumannBC
        variable = p
        boundary = 'Block2_domain'
        v = across_flux_sec
    [../]
    [./pressure_top]
        type = DirichletBC
        variable = p
        boundary = 'top'
        value = 0.0
    [../]
    [./pressure_right]
        type = DirichletBC
        variable = p
        boundary = 'right'
        value = 0.0
    [../]
    [./flux_bot]
        type = NeumannBC
        variable = p
        boundary = '1'
        value = 0.0
    [../]
    [./flux_left]
        type = NeumannBC
        variable = p
        boundary = 'left'
        value = 0.0
    [../]
    [./injection]
        type = NeumannBC
        variable = p
        boundary = '6'
        value = 0.17e-11
    [../]
[]



[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
    petsc_options_value = '200 hypre boomeramg'
  [../]
[]

[Executioner]
    type = Transient
    dt = 0.00004
    automatic_scaling = true
    end_time = 1.0
    [TimeIntegrator]
        type = CentralDifference
    []
[]


[Outputs]
    exodus = true
    interval = 1
[]


[MultiApps]
    #allocate transfer from mainApp to subApp
    [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'sub.i'
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    #get displacement residuals from subApp to mainApp
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'disp_plusminus_sub_scaled_x disp_plusminus_sub_scaled_y traction_sub_strike traction_sub_normal sliprate_sub_strike sliprate_sub_normal slip_sub_strike slip_sub_normal statevar_sub'
        variable = 'disp_plusminus_scaled_x disp_plusminus_scaled_y traction_strike traction_normal sliprate_strike sliprate_normal slip_strike slip_normal statevar'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #push system residual vector from mainApp to subApp
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'resid_x resid_y resid_visco_x resid_visco_y resid_pressure_x resid_pressure_y resid_damp_x resid_damp_y resid_pressdamp_x resid_pressdamp_y fluid_disp_plusminus_x fluid_disp_plusminus_y fluid_vel_plusminus_x fluid_vel_plusminus_y p'
        variable = 'resid_sub_x resid_sub_y resid_visco_sub_x resid_visco_sub_y resid_pressure_sub_x resid_pressure_sub_y resid_damp_sub_x resid_damp_sub_y resid_pressdamp_sub_x resid_pressdamp_sub_y fluid_disp_plusminus_sub_x fluid_disp_plusminus_sub_y fluid_vel_plusminus_sub_x fluid_vel_plusminus_sub_y p_sub_main'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
[]

# [Debug]
#     show_execution_order = ALWAYS
# []
