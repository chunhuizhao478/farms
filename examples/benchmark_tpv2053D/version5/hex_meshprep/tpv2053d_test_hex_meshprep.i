##########################################################
# TPV205 benchmark
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

[Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -22500
      xmax = 22500
      ymin = -22000
      ymax = 0
      zmin = -14000
      zmax = 14000
      nx = 225
      ny = 110
      nz = 140
      subdomain_ids = 1
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z < 0'
      block_id = 2
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z > 0'
      block_id = 3
    []       
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_2
        split_interface = true
        block_pairs = '2 3'
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
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z' 
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all]
        strain = SMALL
        add_variables = true
      [../]
    [../]
  [../]
[]

[Materials]
    [stress]
        type = ComputeLinearElasticStress
    []
    [elasticity]
      type = ComputeIsotropicElasticityTensor
      shear_modulus = 32.04e9
      lambda = 32.04e9
    []
[]

[Executioner]
    type = Steady
[]

[Outputs]
    exodus = true
[]
