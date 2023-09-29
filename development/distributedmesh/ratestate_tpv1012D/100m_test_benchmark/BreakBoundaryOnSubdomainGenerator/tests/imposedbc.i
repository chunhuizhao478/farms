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
    [break_boundary]
        type = BreakBoundaryOnSubdomainGenerator
        input = interface
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
    [disp_plus_scaled_x]
        order = FIRST
        family = LAGRANGE
        initial_condition = 1
    []
    [disp_minus_scaled_x]
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

[BCs]
    [./matchval_primary_x]
        type = MatchedValueBC
        variable = disp_x
        v = disp_plus_scaled_x
        boundary = 'Block0_Block1'
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
