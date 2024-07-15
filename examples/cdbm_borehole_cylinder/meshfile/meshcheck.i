[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  'cylinder_w_hole.msh'
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
    exodus = false
    nemesis = true
[]