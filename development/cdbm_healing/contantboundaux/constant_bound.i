[Mesh]
    type = GeneratedMesh
  
    dim = 2
  
    xmin = 0
    xmax = 1
  
    ymin = 0
    ymax = 1
  
    nx = 10
    ny = 10
  []
  
  [Variables]
    [u]
      order = FIRST
      family = MONOMIAL
    []
  
    [v]
        order = FIRST
        family = MONOMIAL
    []
  []
  
  [AuxVariables]
    [bounds_dummy]
        order = FIRST
        family = MONOMIAL
    []
  []
  
  [Kernels]
    [diff_u]
      type = Diffusion
      variable = u
    []
  
    [diff_v]
      type = Diffusion
      variable = v
    []
  []

  
  [Bounds]
    [u_upper_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = u
      bound_type = upper
      bound_value = 1
    []
    [u_lower_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = u
      bound_type = lower
      bound_value = 0
    []
  
    [v_upper_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = v
      bound_type = upper
      bound_value = 3
    []
    [v_lower_bound]
      type = ConstantBoundsAux
      variable = bounds_dummy
      bounded_variable = v
      bound_type = lower
      bound_value = -1
    []
  []
  
  [Executioner]
    type = Steady
  
    solve_type = 'PJFNK'
    petsc_options_iname = '-snes_type'
    petsc_options_value = 'vinewtonrsls'
  []
  
  [Outputs]
    exodus = true
  []