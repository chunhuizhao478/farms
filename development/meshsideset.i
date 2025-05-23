[Mesh]
construct_node_list_from_side_list = true
    [msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 10
        ny = 10
        xmin = -50
        xmax = 50
        ymin = -50
        ymax = 50
        elem_type = QUAD9
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y > 0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
        add_interface_on_two_sides = true
    []
    [sideset1]
        type = ParsedGenerateSideset
        input = split
        combinatorial_geometry = ' y>9 & y<11'
        new_sideset_name = upper_interior
    []
    [sideset2]
        type = ParsedGenerateSideset
        input = sideset1
        combinatorial_geometry = 'y<-9 & y>-11'
        new_sideset_name = lower_interior
    []
    [add_nodesets2]
        type = NodeSetsFromSideSetsGenerator
        input = sideset2
    []

[]

[GlobalParams]
    displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./nearest_node_value_p]
  [../]
  [./nearest_node_value_m]
  [../]
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            add_variables = true
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

[AuxKernels]
  [./nearest_node_value_p]
    type = NearestNodeValueAux
    variable = nearest_node_value_p
    boundary = 'Block0_Block1'
    paired_variable = disp_x
    paired_boundary = upper_interior
  [../]
  [./nearest_node_value_m]
    type = NearestNodeValueAux
    variable = nearest_node_value_m
    boundary = 'Block1_Block0'
    paired_variable = disp_x
    paired_boundary = lower_interior
  [../]
[]

[BCs]
  [./y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./x_top]
    type = FunctionDirichletBC
    boundary = bottom
    variable = disp_x
    function = disp
  [../]
  [./x_bot]
    type = FunctionDirichletBC
    boundary = top
    variable = disp_x
    function = disp
  [../]
  [./x_bot1]
    type = MatchedValueBC
    boundary = 'Block0_Block1'
    variable = disp_x
    v = nearest_node_value_p
  [../]
  [./x_bot2]
    type = MatchedValueBC
    boundary = 'Block1_Block0'
    variable = disp_x
    v = nearest_node_value_m
  [../]
[]

[Functions]
  [./disp]
    type = PiecewiseLinear
    x = '0.0 1.0 2.0 3.0 4.0' # time
    y = '0.0 1.0 0.0 -1.0 0.0'  # displacement
  [../]
[]

[Materials]
  [./elasticity_tensor_block]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.25
  [../]
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Executioner]
    type = Transient
    dt = 0.005
    end_time = 5
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]


[Outputs]
  exodus = true
[]


