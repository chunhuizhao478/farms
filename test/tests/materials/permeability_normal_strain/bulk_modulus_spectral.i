# Unit test: bulk modulus extracted from the degraded SPECTRAL elastic tangent
# of NDSmallDeformationIsotropicElasticity via K = (1/9) I:C:I.
#
# Geometry: single QUAD4 of size 1 x 1 m on [0,1]^2.
# Strain is imposed kinematically through the disp_x/disp_y AuxVariables
# (set by FunctionAux), so ComputeSmallStrain sees a uniform strain and the
# u-diffusion solve only exists to make the Steady executioner non-empty.
#
# Material constants (E = 50e9, nu = 0.3):
#   K = E/(3(1-2nu)) = 4.1666666667e10
#   G = E/(2(1+nu)) = 1.9230769231e10
# AT1 degradation g(d) = (1-d)^2 (1-eta) + eta, eta = 1e-6.
#
# Closed form (see EXPECTED_VALUES.md):
#   K_eff = K + (g-1) * ( lambda*H(tr eps) + (2G/9)*N+ ),  lambda = K - 2G/3
# with H(tr eps) = 1 if tr eps > 0 else 0 and N+ = #{positive principal strains}.
#
# The d/disp functions and the output file_base are overridden per test case via
# cli_args in the `tests` spec (undamaged / compression / tension).
# Defaults below reproduce the undamaged-tension case (g=1 -> K_eff = K).

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
  # Trivial variable so the Steady solve is non-empty.
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
    order = FIRST
    family = LAGRANGE
  []
  [Kd]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [disp_x_aux]
    type = FunctionAux
    variable = disp_x
    function = '1e-3 * x'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [disp_y_aux]
    type = FunctionAux
    variable = disp_y
    function = '1e-3 * y'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Kd]
    type = MaterialRealAux
    variable = Kd
    property = bulk_modulus_degraded
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
  [u_left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [u_right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 50e9
    poissons_ratio = 0.3
  []
  [strain]
    type = ComputeSmallStrain
  []
  [bulk]
    type = GenericConstantMaterial
    prop_names = 'K G'
    prop_values = '4.1666666667e10 1.9230769231e10'
  []
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    decomposition = SPECTRAL
    model_type = AT1
    eta = 1e-6
    porous_flow_coupling = false
  []
  [stress]
    type = NDComputeSmallDeformationStress
    elasticity_model = elasticity
  []
[]

[Postprocessors]
  [Kd_avg]
    # Unit element => element average == the single QP value.
    type = ElementAverageValue
    variable = Kd
    execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
    # Explicit file_base so the produced file is exactly <file_base>.csv
    # (overridden per test case via cli_args).
    file_base = bulk_modulus_spectral
  []
[]
