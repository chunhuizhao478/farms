##########################################################
# TPV205 benchmark
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

##########################################################
# TPV205 benchmark
# Global coordinate is the same as local coordinate
# x - strike; y - dip; z - normal
##########################################################

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -2
    xmax = 2
    ymin = -4
    ymax = 0
    zmin = -2
    zmax = 2
    nx = 2
    ny = 2
    nz = 2
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= -1 & x <= 1 & y >= -4 & z < 0 '
    block_id = 5
  []
  [./new_block_2]
    type = ParsedSubdomainMeshGenerator
    input = new_block_1
    combinatorial_geometry = 'x >= -1 & x <= 1 & y >= -4 & z > 0'
    block_id = 6
  []       
  [./split_1]
    type = BreakMeshByBlockGenerator
    input = new_block_2
    split_interface = true
    block_pairs = '5 6'
  []      
[]

[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[Kernels]
    [./inertia_x]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_x
    []
    [./inertia_y]
      type = NullKernel
      use_displaced_mesh = false
      variable = disp_y
    []
    [./inertia_z]
      type = NullKernel
      use_displaced_mesh = false
      variable = disp_z
    []
[]

[Materials]
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 1
    []
[]

[Executioner]
    type = Transient
    dt = 1
    end_time = 12.0
    num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
[]