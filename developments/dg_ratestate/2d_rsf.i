#Apr27, 2024

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 80
        ny = 8
        xmin = -10
        xmax = 10
        ymin = -1
        ymax = 1
        elem_type = QUAD4
    [../]
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    [../]
    [./mark_fault]
        type = SideSetsBetweenSubdomainsGenerator
        input = new_block
        new_boundary = fault
        primary_block = 0
        paired_block = 1
    [../]
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
[]

[Physics/SolidMechanics/QuasiStatic]
    [all]
      add_variables = true
    []
[]

[Materials]
    [elasticity]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 1e9
      poissons_ratio = 0.3
    []
    [stress]
      type = ComputeLinearElasticStress
    []
[]

[Executioner]
    type = Transient
    end_time = 5
    dt = 1
[]
  
[Outputs]
    exodus = true
[]

