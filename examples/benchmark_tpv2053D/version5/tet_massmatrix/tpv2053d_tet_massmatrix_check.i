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
      xmin = -4
      xmax = 4
      ymin = -8
      ymax = 0
      zmin = -4
      zmax = 4
      nx = 4
      ny = 4
      nz = 4
      elem_type = TET4
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'z < 0 '
      block_id = 5
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'z > 0'
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
        type = NullKernel
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

  [UserObjects]
    [reader_node]
        type = PropertyReadFile
        prop_file_name = 'massbyidpostprocess.txt'
        read_type = 'node'
        nprop = 1 # number of columns in CSV
    []
    [./nodal_area]
        type = NodalArea
        variable = nodearea
        boundary = 'Block5_Block6'
        execute_on = 'initial timestep_end'
    [../]
  []

  [Functions]
    [node]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'reader_node'
        read_type = 'node'
        # 0-based indexing
        column_number = '0'
    []    
  []

  [AuxVariables]
    [nodemass]
        order = FIRST
        family = LAGRANGE
    []
    [nodearea]
        order = FIRST
        family = LAGRANGE        
    []
  []

  [AuxKernels]
    [nodemass]
        type = FunctionAux
        variable = nodemass
        function = node
    []
  []
  
  [Outputs]
      exodus = true
  []