# Error-path test: the STRAIN porosity law requires the kinematic strain tensor
# (mechanical_strain). This input selects porosity_update_model = strain but provides
# NO ComputeSmallStrain material, so binding `strain_property` must fail loudly rather
# than silently using a zero strain. Paired with a RunException entry expecting an error
# mentioning 'mechanical_strain'.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [d]
    order = CONSTANT
    family = MONOMIAL
  []
  # Consumes the damaged-porosity property so the material is active and its
  # mechanical_strain request is enforced (otherwise MOOSE prunes the unused
  # material and the missing-strain error never fires).
  [porosity_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [porosity]
    type = MaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'INITIAL'
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Materials]
  # NOTE: no ComputeSmallStrain here on purpose -> mechanical_strain is undefined.
  [porosity_strain]
    type = ElkPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_update_model = strain
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
