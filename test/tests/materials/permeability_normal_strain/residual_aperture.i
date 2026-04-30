# Unit test: residual aperture w_r (Heider 2021 eq. 46 closed-crack branch).
#
# Verifies the closed-branch term of
#   w_h = max{ (f_c * w_c) * chi_d,   (open)
#              (f_c * w_r) * chi_d }  (closed)
# When the open-crack aperture w_c is much smaller than w_r, the `max` selects
# the closed branch and K_frac = (f_c*w_r)^2/12 * (I - n_d⊗n_d) instead of the
# tiny w_c-driven value.
#
# Setup:
#   - disp_x = 1e-3 * x, disp_y = 0  ->  eps_xx = 1e-3 (unique most-tensile
#     eigenvalue of the strain tensor), all other components zero. The
#     principal-strain fallback therefore returns n_d = e_x deterministically
#     on every platform (not dependent on LAPACK tie-breaking of a degenerate
#     zero tensor).
#   - characteristic_length_value = h_c = 1e-7 m, so the open-crack aperture
#     w_c = h_c * (1 + eps_nn) ~ 1.001e-7 m << w_r = 1e-5 m, forcing the
#     residual (closed) branch of Heider eq. (46) to win.
#   - Uniform d = 0.7 -> grad(d) = 0 -> crack_normal_source = principal_strain.
#   - residual_aperture = 1e-5 m, fc = 1, anisotropic = true.
#
# Hand calc:
#   eps_nn = e_x . eps . e_x = eps_xx = 1e-3
#   w_c = h_c * (1 + 1e-3) = 1.001e-7
#   w_h = max(1*1.001e-7*1, 1*1e-5*1) = 1e-5  <- closed branch wins
#   k_w = (1e-5)^2 / 12 = 8.3333e-12
#   d^b = 0.7^2 = 0.49
#   K_frac = k_w * (I - e_x ⊗ e_x) = diag(0, k_w, k_w)
#   K = k0*I + 0.49 * K_frac
#     K_xx        = k0            = 5e-19
#     K_yy = K_zz = k0 + 0.49*k_w = 4.0833e-12
#     K_xy        = 0
#
# If residual_aperture defaulted to 0 (legacy behavior), w_h would reduce to
# f_c*w_c*chi_d ~ 1.001e-7 m, k_w = (1.001e-7)^2/12 ~ 8.35e-16, and
# K_yy = K_zz would be ~4.09e-16 -- about 4 orders of magnitude smaller than
# the residual-branch value computed above. The Exodiff against the new gold
# therefore distinguishes residual_aperture=1e-5 from residual_aperture=0
# by ~10^4.

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
    function = '1e-3 * x'   # eps_xx = 1e-3, unique most-tensile eigenvalue
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
    crack_normal_source = principal_strain  # grad(d) = 0 for uniform d
    characteristic_length_type = constant
    characteristic_length_value = 1e-7      # smaller than w_r to force closed branch
    permeability_anisotropic = true
    damage_threshold_for_permeability = 0.5
    correction_factor_fc = 1.0
    residual_aperture = 1e-5                # closed-branch floor (Heider eq. 46)
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
