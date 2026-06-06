# Unit test: Heider (2021) normal-strain permeability with the REGULARIZED
# crack normal (regularize_crack_normal = true).
#
# Purpose: verify that at the fully-damaged crack core, where grad(d) = 0, the
# regularized normal n_d = grad(d) / (|grad(d)| + eps) collapses to 0, so the
# tangential projector (I - n_d (x) n_d) -> I and the fracture permeability is
# ISOTROPIC there -- instead of the legacy hard-cutoff behavior, which would
# detect |grad(d)| < tol and fall back to the matrix permeability K = k0*I.
#
# Geometry: single QUAD4 of size 1 x 1 m on [0,1]^2.
# Damage:   uniform d = 0.7  ->  grad(d) = 0 exactly at every QP.
# Strain:   disp_x = 1e-3 * x, disp_y = 0  ->  eps_xx = 1e-3, others 0.
#           (Irrelevant to the result: with n_d = 0, eps_nn = n_d.eps.n_d = 0.)
#
# Because grad(d) = 0, n_d = (0,0,0)/(0 + eps) = (0,0,0) for ANY eps > 0:
#   eps_nn = 0
#   w_c    = h_c * |1 + eps_nn| = h_c = 1e-3
#   chi_d  = H(d - 0.5) = 1            (d = 0.7 >= threshold)
#   w_h    = f_c * w_c * chi_d = 1e-3
#   k_w    = w_h^2 / 12 = (1e-3)^2 / 12 = 8.333333e-8
#   (I - n_d (x) n_d) = I              (n_d = 0)  -> ISOTROPIC even though
#                                       permeability_anisotropic = true
#   d^b    = 0.7^2 = 0.49              (uniform, no QP averaging)
# Expected (isotropic, despite permeability_anisotropic = true):
#   K_xx = K_yy = K_zz = k0 + 0.49 * k_w = 5e-19 + 0.49 * 8.333333e-8
#                      = 4.083333e-8
#   K_xy = 0
#
# Contrast with the legacy (regularize_crack_normal = false) path on this same
# setup: |grad(d)| = 0 < damage_gradient_tolerance -> have_normal = false ->
# K = k0*I = 5e-19 (matrix perm). The ~11-order-of-magnitude jump between the
# two paths is what this test pins down.

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
    order = FIRST
    family = LAGRANGE
  []
  [eff_perm_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_22]
    order = CONSTANT
    family = MONOMIAL
  []
  [eff_perm_01]
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
    function = '0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_aux]
    type = FunctionAux
    variable = d
    function = '0.7'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [eff_perm_00]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 0
    variable = eff_perm_00
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_11]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 1
    column = 1
    variable = eff_perm_11
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_22]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 2
    column = 2
    variable = eff_perm_22
    execute_on = 'TIMESTEP_END'
  []
  [eff_perm_01]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 1
    variable = eff_perm_01
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
    porous_flow_coupling = true
    permeability_model = normal_strain
    intrinsic_permeability = 5e-19
    perm_exponent = 2
    crack_normal_source = damage_gradient
    # New regularized-normal option under test:
    regularize_crack_normal = true
    crack_normal_regularization = 1e-6   # value is immaterial here: n_d = 0/(0+eps) = 0
    characteristic_length_type = constant
    characteristic_length_value = 1e-3
    permeability_anisotropic = true      # yet the result is isotropic because n_d = 0
    damage_threshold_for_permeability = 0.5
    correction_factor_fc = 1.0
  []
  [stress]
    type = NDComputeSmallDeformationStress
    elasticity_model = elasticity
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
[]
