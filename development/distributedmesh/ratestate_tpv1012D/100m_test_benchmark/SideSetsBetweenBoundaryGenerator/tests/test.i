[Mesh]
    [generated0]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 201
        ny = 200
        xmin = -10050
        xmax = 10050
        ymin = 0
        ymax = 10000
        boundary_name_prefix = pillar1
    []
    [rename_block0]
        type = RenameBoundaryGenerator
        input = generated0
        old_boundary = '0 1 2 3'
        new_boundary = 'block0_bottom block0_right block0_left block0_top'
    []
    [block0_id]
        type = SubdomainIDGenerator
        input = rename_block0
        subdomain_id = 0
    []
    [generated1]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 201
        ny = 200
        xmin = -10050
        xmax = 10050
        ymin = -10000
        ymax = 0
        boundary_name_prefix = pillar2
        boundary_id_offset = 4
    []
    [rename_block1]
        type = RenameBoundaryGenerator
        input = generated1
        old_boundary = '4 5 6 7'
        new_boundary = 'block1_bottom block1_right block1_left block1_top'
    []
    [block1_id]
        type = SubdomainIDGenerator
        input = rename_block1
        subdomain_id = 1
    []
    [collect]
        type = MeshCollectionGenerator
        inputs = 'block0_id block1_id'
    []
    [sideset]
        type = SideSetsBetweenSubdomainsGenerator
        input = collect
        primary_block = 0
        paired_block = 1
        new_boundary = 'interface'
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'

[]

[Variables]
    [./disp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    #interface displacement boundary condition (scale disp_plusminus_(x/y))
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

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_x disp_y'
            planar_formulation = PLANE_STRAIN
        [../]
        [../]
    [../]
[]

[Kernels]
    [./inertia_x]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_x
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y
    []
[]

[InterfaceKernels]
    #apply displacement prediction and retrieve its residuals
    [./ratestate_x]
        type = RateStateInterfaceKernelGlobalxdev
        variable = disp_x
        neighbor_var = disp_x
        disp_strike_plus = disp_plus_scaled_x
        disp_strike_minus = disp_minus_scaled_x
        boundary = 'interface'
    []
    [./ratestate_y]
        type = RateStateInterfaceKernelGlobalydev
        variable = disp_y
        neighbor_var = disp_y
        disp_normal_plus = disp_plus_scaled_y
        disp_normal_minus = disp_minus_scaled_y
        boundary = 'interface'
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
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
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