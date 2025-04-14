# Mesh Check #

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = 'tpv2053d_200m.msh'
    []
    [./new_block_1]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'x >= -15000 & x <= 15000 & z >= -15000 & y < 0'
        block_id = 100
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x >= -15000 & x <= 15000 & z >= -15000 & y > 0'
        block_id = 200
    []       
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_2
        split_interface = true
        block_pairs = '100 200'
    []      
    [./sidesets]
      input = split_1
      type = SideSetsFromNormalsGenerator
      normals = '-1 0 0
                  1 0 0
                  0 -1 0
                  0 1 0
                  0 0 -1
                  0 0 1'
      new_boundary = 'left right bottom top back front'
    [] 
[]

[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
[]

[Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          add_variables = true
        []
      []
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
[]

[Executioner]
    type = Steady
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
    nl_rel_tol = 1e-12
    nl_abs_tol = 1e-50
[]

[Outputs]
    exodus = true
[]

[AuxVariables]
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
[]

# [UserObjects]
#     [./nodal_area]
#         type = NodalArea
#         variable = nodal_area
#         boundary = 'interface'
#         execute_on = 'initial TIMESTEP_BEGIN'
#     [../]
# []