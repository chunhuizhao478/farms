[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 801
        ny = 800
        xmin = -10012.5
        xmax = 10012.5
        ymin = -10000
        ymax = 10000
    [../]
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    #add "Block0_Block1" and "Block1_Block0" interfaces
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
        add_interface_on_two_sides = true
    []
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y'
    
    ##element length (m)
    len = 25
    
    ##rate-and-state coefficients
    f_o = 0.6
    rsf_a = 0.008
    rsf_b = 0.012
    rsf_L = 0.02
    delta_o = 1e-6

    ##initial normal traction (Pa)
    Tn_o = 120e6

    ##initial shear traction (Pa)
    Ts_o = 75e6

    ##initial sliprate (m/s)
    Vini = 1e-12

    ##initial state variable
    statevarini = 1.606238999213454e9
[]

[Variables]
    [disp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    #fault displacement residual from interfacekernel
    [disp_plusminus_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
    [disp_plusminus_sub_scaled_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_plusminus_sub_scaled_y]
        order = FIRST
        family = LAGRANGE
    []
    #side element volume 
    [element_side_volume]
        order = FIRST
        family = LAGRANGE
    []
    #residual received from mainApp
    [resid_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [resid_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    #residual received from mainApp (damping)
    [resid_damp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [resid_damp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    ##initial shear stress
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []
    #traction
    [traction_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #sliprate
    [sliprate_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [sliprate_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #slip
    [slip_sub_strike]
        order = CONSTANT
        family = MONOMIAL
    []
    [slip_sub_normal]
        order = CONSTANT
        family = MONOMIAL
    []
    #statevar
    [statevar_sub]
        order = CONSTANT
        family = MONOMIAL
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
        value = 25
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
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_traction_normal]
        type = MaterialRealAux
        property = traction_normal
        variable = traction_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_sliprate_strike]
        type = MaterialRealAux
        property = sliprate_strike
        variable = sliprate_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_sliprate_normal]
        type = MaterialRealAux
        property = sliprate_normal
        variable = sliprate_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_strike]
        type = MaterialRealAux
        property = slip_strike
        variable = slip_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip_normal]
        type = MaterialRealAux
        property = slip_normal
        variable = slip_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_statevar]
        type = MaterialRealAux
        property = statevar
        variable = statevar_sub
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_sub_x disp_sub_y'
            planar_formulation = PLANE_STRAIN
        [../]
        [../]
    [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
    [../]
[]

[InterfaceKernels]
    #apply displacement prediction and retrieve its residuals
    [./ratestate_x]
        type = RateStateInterfaceKernelGlobalx
        variable = disp_sub_x
        neighbor_var = disp_sub_x
        extra_vector_tags = 'restore_tag_x'
        boundary = 'Block0_Block1'
    []
    [./ratestate_y]
        type = RateStateInterfaceKernelGlobaly
        variable = disp_sub_y
        neighbor_var = disp_sub_y
        extra_vector_tags = 'restore_tag_y'
        boundary = 'Block0_Block1'
    []
[]

[Problem]
    extra_tag_vectors = 'restore_tag_x restore_tag_y'
[]

[Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    [./czm_mat]
        type = RateStateFrictionLaw2DAsBC
        reaction_rsf_x  = resid_sub_x
        reaction_rsf_y  = resid_sub_y
        reaction_damp_x = resid_damp_sub_x
        reaction_damp_y = resid_damp_sub_y
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'Block0_Block1'
        output_properties = 'sliprate_strike slip_strike statevar traction_strike traction_normal alongfaultdisp_strike_plus alongfaultdisp_strike_minus'
        # outputs = exodus
    [../]
[]

[UserObjects]
    #compute element side volume (using CONTACT modulus)
    # [element_side_volume]
    #     type = NodalArea
    #     variable = element_side_volume
    #     boundary = 'Block0_Block1 Block1_Block0'
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

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 20
[]