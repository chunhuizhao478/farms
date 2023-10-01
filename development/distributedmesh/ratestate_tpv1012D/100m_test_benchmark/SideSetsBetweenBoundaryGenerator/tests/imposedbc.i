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
    [interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = new_block
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

[Variables]
    [./disp_x_top]
        order = FIRST
        family = LAGRANGE
        block = 0
    [../]
    [./disp_x_bottom]
        order = FIRST
        family = LAGRANGE
        block = 1
    [../]
    [./disp_y_top]
        order = FIRST
        family = LAGRANGE
        block = 0
    [../]
    [./disp_y_bottom]
        order = FIRST
        family = LAGRANGE
        block = 1
    [../]
[]

[AuxVariables]
    [disp_plus_scaled_x]
        order = FIRST
        family = LAGRANGE
        initial_condition = 1
    []
    [disp_minus_scaled_x]
        order = FIRST
        family = LAGRANGE
        initial_condition = -1
    []
    [disp_plus_scaled_y]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
    []
    [disp_minus_scaled_y]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
    []
[]

[Kernels]
    [./disp_top_x]
        type = StressDivergenceTensors
        variable = disp_x_top
        displacements = 'disp_x_top disp_y_top'
        component = 0
    []
    [./disp_bot_x]
        type = StressDivergenceTensors
        variable = disp_x_bottom
        displacements = 'disp_x_bottom disp_y_bottom'
        component = 0
    []
    [./disp_top_y]
        type = StressDivergenceTensors
        variable = disp_y_top
        displacements = 'disp_x_top disp_y_top'
        component = 1
    []
    [./disp_bot_y]
        type = StressDivergenceTensors
        variable = disp_y_bottom
        displacements = 'disp_x_bottom disp_y_bottom'
        component = 1
    []
    [./inertia_x_top]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x_top
        block = 0
    []
    [./inertia_y_top]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y_top
        block = 0
    []
    [./inertia_x_bot]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x_bottom
        block = 1
    []
    [./inertia_y_bot]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y_bottom
        block = 1
    []
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
    [strain_top]
        type = ComputeSmallStrain
        displacements = 'disp_x_top disp_y_top'
        block = 0
    []
    [strain_bottom]
        type = ComputeSmallStrain
        displacements = 'disp_x_bottom disp_y_bottom'
        block = 1
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
[]


# [InterfaceKernels]
#     #apply displacement prediction and retrieve its residuals
#     [./ratestate_x]
#         type = RateStateInterfaceKernelGlobalxdev
#         variable = disp_x_top
#         neighbor_var = disp_x_bottom
#         disp_strike_plus = disp_plus_scaled_x
#         disp_strike_minus = disp_minus_scaled_x
#         boundary = 'Block0_Block1'
#     []
#     [./ratestate_y]
#         type = RateStateInterfaceKernelGlobalydev
#         variable = disp_y_top
#         neighbor_var = disp_y_bottom
#         disp_normal_plus = disp_plus_scaled_y
#         disp_normal_minus = disp_minus_scaled_y
#         boundary = 'Block0_Block1'
#     []
# []

[BCs]
    #assign displacement boundary condition
    [./matchval_primary_xt]
        type = MatchedValueBC
        variable = disp_x_top
        v = disp_plus_scaled_x
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_xb]
        type = MatchedValueBC
        variable = disp_x_bottom
        v = disp_minus_scaled_x
        boundary = 'Block1_Block0'
    []
    [./matchval_primary_yt]
        type = MatchedValueBC
        variable = disp_y_top
        v = disp_plus_scaled_y
        boundary = 'Block0_Block1'
    []
    [./matchval_secondary_yb]
        type = MatchedValueBC
        variable = disp_y_bottom
        v = disp_minus_scaled_y
        boundary = 'Block1_Block0'
    []
[]

[Executioner]
    type = Transient
    dt = 1
    end_time = 5.0
    num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]
