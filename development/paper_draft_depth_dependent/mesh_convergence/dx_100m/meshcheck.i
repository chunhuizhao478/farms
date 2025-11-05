##main fault parameters
xmin_fault = -15000 #xmin of fault
xmax_fault = 15000 #xmax of fault
ymin_domain = -2000 #ymin of domain
ymax_domain = 2000 #ymax of domain

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/2dmesh_100m.msh'
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y>0 & x>${xmin_fault} & x<${xmax_fault} & y<${ymax_domain}'
      block_id = 100
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'y<0 & x>${xmin_fault} & x<${xmax_fault} & y>${ymin_domain}'
      block_id = 200
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block_2
      split_interface = true
      block_pairs = '100 200'
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

[AuxVariables]
    [./ini_shear_sts]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_normal_sts]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_d]
        order = CONSTANT
        family = MONOMIAL
    []
    [./D_c]
        order = CONSTANT
        family = MONOMIAL
    []
[]
