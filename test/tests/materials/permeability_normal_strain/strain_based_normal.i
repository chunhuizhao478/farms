# Unit test: strain-based crack normal (Liu et al. 2024 CMAME, eqs. 29-30).
#
# Purpose: verify that with crack_normal_source = principal_strain the crack
# normal is n_F = e_1, the eigenvector of the MAXIMUM principal strain of the
# model's mechanical strain -- NOT the damage gradient. The test is constructed
# so the two normals are ORTHOGONAL, so a strain-based normal and a
# damage-gradient normal give different effective-permeability tensors.
#
# Geometry: single QUAD4 of size 1 x 1 m on [0,1]^2 (2x2 Gauss).
# Damage:   d(x,y) = 0.5 + 0.4*x  ->  grad(d) = (0.4, 0, 0)  (|| e_x).
#           A damage-gradient normal would therefore be n_d = e_x.
# Strain:   disp_x = 0, disp_y = 1e-3*y  ->  eps = diag(0, 1e-3, 0).
#           Max principal strain eps_1 = 1e-3 with eigenvector e_1 = e_y.
#           So the STRAIN-BASED normal is n_F = e_y, orthogonal to grad(d).
#
# Hand calc (strain-based normal n_F = e_y, Liu eqs. 29-30):
#   eps_nn = n_F . eps . n_F = eps_yy = 1e-3
#   w_c    = h_c * |1 + eps_nn| = 1e-3 * 1.001 = 1.001e-3
#   chi_d  = H(d - 0.5) = 1                 (d >= 0.5 at every QP)
#   w_h    = f_c * w_c * chi_d = 1.001e-3
#   k_w    = w_h^2 / 12 = (1.001e-3)^2 / 12 = 8.350008e-8
#   (I - n_F (x) n_F) = diag(1, 0, 1)       (n_F = e_y) -> blocks y, enhances x,z
#   d_qp   = 0.5 + 0.4*{0.2113, 0.7887} = {0.58452, 0.81548}
#   <d^2>  = 1/2*(0.58452^2 + 0.81548^2) = 0.503336
#   alpha  = <d^2> * k_w = 0.503336 * 8.350008e-8 = 4.202838e-8
# Expected:
#   K_yy = k0 = 5e-19            (NORMAL direction blocked -> the discriminator)
#   K_xx = K_zz = k0 + alpha = 4.202838e-8
#   K_xy = 0
#
# Discriminator: a damage-gradient normal (n_d = e_x) would instead block x,
# giving K_xx = k0 and K_yy = K_zz = alpha. The swap of which diagonal entry
# collapses to k0 proves the normal is strain-based (e_1), not damage-based.

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
    function = '0'
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
    function = '0.5 + 0.4 * x'
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
    # Strain-based crack normal under test: n_F = e_1 (Liu 2024 eqs. 29-30).
    crack_normal_source = principal_strain
    characteristic_length_type = constant
    characteristic_length_value = 1e-3
    permeability_anisotropic = true
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
