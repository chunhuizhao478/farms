[Mesh]
    [./msh]
        type = FileMeshGenerator
        # file =  './planarfault3D.msh'
        file = '../TPV24.msh'
    []
    [./new_block_1]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'x > -16000 & x < 12000 & y < 0 & z > -15000'
        block_id = 1
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x > -16000 & x < 12000 & y > 0 & z > -15000'
        block_id = 2
    []    
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_2
        split_interface = true
        block_pairs = '1 2'
    []
    [./new_block_3]
        type = ParsedSubdomainMeshGenerator
        input = split_1
        combinatorial_geometry = '0.5773505384 * x + y < 0 & z > -15000 & x >= -100 & y <= 100 & x <= 10500 & y >= -6500'
        block_id = 3
    []
    [./new_block_4]
        type = ParsedSubdomainMeshGenerator
        input = new_block_3
        combinatorial_geometry = '0.5773505384 * x + y > 0 & z > -15000 & x >= -100 & y <= 100 & x <= 10500 & y >= -6500'
        block_id = 4
    []    
    [./split_2]
        type = BreakMeshByBlockGenerator
        input = new_block_4
        split_interface = true
        block_pairs = '3 4'
    []   
[]

[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            add_variables = true
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