##########################################################
# TPV205 benchmark
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

[Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -32000
      xmax = 32000
      ymin = -64000
      ymax = 0
      zmin = -32000
      zmax = 32000
      nx = 20
      ny = 20
      nz = 20
      subdomain_ids = 1
    []
    [./fault_area_block_1]
      type = SubdomainBoundingBoxGenerator
      input = msh
      block_id = 2
      bottom_left = '-16000 -16000 3200'
      top_right = '16000 0 -3200'
      location = INSIDE
    []
    [./refine_fault_area_block_1]
      type = RefineBlockGenerator
      input = fault_area_block_1
      block = '2'
      refinement = '1'
      enable_neighbor_refinement = false
    []
    [./fault_area_block_2]
      type = SubdomainBoundingBoxGenerator
      input = refine_fault_area_block_1
      block_id = 3
      bottom_left = '-15200 -15200 800'
      top_right = '15200 0 -800'
      location = INSIDE
    []
    [./refine_fault_area_block_2]
      type = RefineBlockGenerator
      input = fault_area_block_2
      block = '3'
      refinement = '1'
      enable_neighbor_refinement = false
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = refine_fault_area_block_2
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z < 0 & z > -800'
      block_id = 4
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y >= -15000 & z > 0 & z < 800'
      block_id = 5
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
