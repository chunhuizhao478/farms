[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../gen/refinedmesh.msh'
    []
    [./subdomain_id] 
        input = msh
        type = SubdomainPerElementGenerator 
        element_ids = 
        
        
        subdomain_ids = 
        
        
    []
    [./split]
        input = subdomain_id
        type = BreakMeshByBlockGenerator
        surrounding_blocks = 
        
        
        
        interface_name = czm
    []
    [./sidesets]
        input = split
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            add_variables = true
            planar_formulation = PLANE_STRAIN
            generate_output = 'stress_xx stress_yy stress_xy'
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