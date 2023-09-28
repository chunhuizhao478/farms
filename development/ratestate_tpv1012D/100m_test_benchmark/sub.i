[Mesh]
    [generated1]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 201
        ny = 200
        xmin = -10050
        xmax = 10050
        ymin = -10000
        ymax = 10000
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = generated1
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    [break_boundary]
        type = BreakBoundaryOnSubdomainGenerator
        input = new_block
    []
    [interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = break_boundary
        primary_block = 0
        paired_block = 1
        new_boundary = 'Block0_Block1'
    []
    [interface2]
        type = SideSetsBetweenSubdomainsGenerator
        input = interface
        primary_block = 1
        paired_block = 0
        new_boundary = 'Block1_Block0'
    []
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y'
    
    ##element length (m)
    len = 100
    
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

    #
    f_w = 0.2
    Vw = 0.1
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
    #alongfaultdisp
    [alongfaultdisp_strike_plus_sub]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_strike_minus_sub]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_normal_plus_sub]
        order = CONSTANT
        family = MONOMIAL
    []
    [alongfaultdisp_normal_minus_sub]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
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
    [output_alongdisp_strike_plus]
        type = MaterialRealAux
        property = alongfaultdisp_strike_plus
        variable = alongfaultdisp_strike_plus_sub
        boundary = 'Block0_Block1 Block1_Block0'
        execute_on = 'TIMESTEP_END'
    []
    [output_alongdisp_strike_minus]
        type = MaterialRealAux
        property = alongfaultdisp_strike_minus
        variable = alongfaultdisp_strike_minus_sub
        boundary = 'Block0_Block1 Block1_Block0'
        execute_on = 'TIMESTEP_END'
    []
    [output_alongdisp_normal_plus]
        type = MaterialRealAux
        property = alongfaultdisp_normal_plus
        variable = alongfaultdisp_normal_plus_sub
        boundary = 'Block0_Block1 Block1_Block0'
        execute_on = 'TIMESTEP_END'
    []
    [output_alongdisp_normal_minus]
        type = MaterialRealAux
        property = alongfaultdisp_normal_minus
        variable = alongfaultdisp_normal_minus_sub
        boundary = 'Block0_Block1 Block1_Block0'
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
        type = IMRateStateFrictionLaw2D
        RSFlaw = 1
        reaction_rsf_x  = resid_sub_x
        reaction_rsf_y  = resid_sub_y
        reaction_damp_x = resid_damp_sub_x
        reaction_damp_y = resid_damp_sub_y
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'Block0_Block1 Block1_Block0'
    [../]
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