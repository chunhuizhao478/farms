#sample test geometry
[Mesh]
    [./msh]
      type = FileMeshGenerator
      file =  '../../meshgenerator/cdbm/contmf/tria/contmfsmall.msh'
    []
    [./new_block]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y<0'
      block_id = 1
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block
      split_interface = true
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
  
  [Variables]
    [alpha_subsub]
      order = FIRST
      family = LAGRANGE
    []
  
    [B_subsub]
      order = FIRST
      family = LAGRANGE
    []
  []
  
  [AuxVariables]
    [bounds_dummy]
      order = FIRST
      family = LAGRANGE
    []
  []
  
  [Kernels]
    [null_alpha]
      type = NullKernel
      variable = alpha_subsub
    []
  
    [null_B]
      type = NullKernel
      variable = B_subsub
    []
  []
  
  [Bounds]
    [alpha_upper_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = alpha_subsub
      bound_type = upper
      bound_value = 1
    []
    [alpha_lower_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = alpha_subsub
      bound_type = lower
      bound_value = 0
    []
  
    [B_upper_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = B_subsub
      bound_type = upper
      bound_value = 1
    []
    [B_lower_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = B_subsub
      bound_type = lower
      bound_value = 0
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = 'PJFNK'
    petsc_options_iname = '-snes_type'
    petsc_options_value = 'vinewtonrsls'
  []

[Outputs]
  exodus = true
[]
