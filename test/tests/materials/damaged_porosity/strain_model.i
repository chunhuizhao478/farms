# Analytic test for ElkPorousFlowDamagedPorosity, STRAIN update law (Liu et al. 2024,
# CMAME 429:117165, eq. 40):  phi(eps) = phi_0 + eps_1,  clamped to [lower, upper],
# where eps_1 is the maximum (most-tensile) principal strain of mechanical_strain.
#
# A prescribed displacement gives a KNOWN, element-uniform strain. disp_x is a
# continuous piecewise-linear field whose slope (= eps_xx) is constant in each of the
# three elements; disp_y = 0 so eps_yy = eps_xy = 0 and (2D plane strain) eps_zz = 0.
# Hence the strain tensor per element is diag(eps_xx, 0, 0) and
#   eps_1 = max(eps_xx, 0, 0).
#
# Element slopes:  elem0 = 1e-3,  elem1 = 2.0,  elem2 = -0.5.
# phi_0 = 0.008, lower = 0.008, upper = 0.999:
#   elem 0: eps_1 = 1e-3  -> phi = 0.009          (eq. 40, in bounds)
#   elem 1: eps_1 = 2.0   -> phi_raw = 2.008      -> clamp upper -> 0.999
#   elem 2: eps_1 = 0     -> phi = 0.008          (compression -> initial porosity)
# See EXPECTED_VALUES.md for the derivation. The strain_lower_clamp test reuses this
# input with porosity_lower_bound = 0.05 to exercise the lower clamp.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 0.1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  # Trivial field so the Steady solve is non-empty.
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
  [porosity_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  # Continuous piecewise-linear displacement; nodal values are the running integral of
  # the per-element slopes (1e-3, 2.0, -0.5) over element length 1/3. Breakpoints align
  # with mesh nodes (x = 1/3, 2/3), so the element gradient equals the slope exactly.
  [disp_x_func]
    type = PiecewiseLinear
    x = '0 0.333333333333 0.666666666667 1'
    y = '0 0.000333333333 0.667000000000 0.500333333333'
  []
[]

[AuxKernels]
  [disp_x_aux]
    type = FunctionAux
    variable = disp_x
    function = disp_x_func
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y_aux]
    type = FunctionAux
    variable = disp_y
    function = '0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [porosity]
    type = MaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'TIMESTEP_END'
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
  # Kinematic strain (declares mechanical_strain consumed by the porosity material).
  [strain]
    type = ComputeSmallStrain
  []
  [porosity_strain]
    type = ElkPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
    porosity_update_model = strain
    # strain_property defaults to 'mechanical_strain'
  []
[]

[Postprocessors]
  [phi_tension]
    type = PointValue
    variable = porosity_aux
    point = '0.1666667 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [phi_clamp_upper]
    type = PointValue
    variable = porosity_aux
    point = '0.5 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [phi_compression]
    type = PointValue
    variable = porosity_aux
    point = '0.8333333 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  execute_on = 'FINAL'
[]
