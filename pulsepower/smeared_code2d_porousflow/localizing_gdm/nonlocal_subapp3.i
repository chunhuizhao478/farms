#initial damage box 1
bottom_left1 = '-0.0025 -2e-4 0'
top_right1 = '0.0025 2e-4 0'

#initial damage box 2
bottom_left2 = '-2e-4 -0.0025 0'
top_right2 = '2e-4 0.0025 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../2dmeshfile/fieldscale_test1_2d_small.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  [./subdomain_id]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left1}
    top_right = ${top_right1}
    location = INSIDE
    block_id = 1
    input = extranodeset1
  []
  [./subdomain_id2]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left2}
    top_right = ${top_right2}
    location = INSIDE
    block_id = 1
    input = subdomain_id
  []
[]

[Variables]
  [nonlocal_eqstrain]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [eqstrain_local]  # Received from main app
    order = CONSTANT
    family = MONOMIAL
  []
  [crack_damage_aux]  # Received from main app
    order = FIRST
    family = MONOMIAL
  []
[]

[Kernels]
  [react_nonlocal]
    type = Reaction
    variable = nonlocal_eqstrain
    rate = 1.0
  []
  [diffusion_nonlocal]
    type = LocalizingCoefDiffusion
    variable = nonlocal_eqstrain
    coef = ${fparse l*l}
    R = 0.02
    eta = 5
  []
  [reaction_local]
    type = CoupledElkFixedLocalEqstrainForce
    variable = nonlocal_eqstrain
    eqstrain_local = eqstrain_local
  []    
[]

[Materials]
  [./crack_damage]
    type = ParsedMaterial
    property_name = crack_damage
    coupled_variables = 'crack_damage_aux'
    expression = 'crack_damage_aux'
    outputs = exodus
  [../]
[]

# [Kernels]
#   [react_nonlocal]
#     type = CoupledReaction
#     variable = nonlocal_eqstrain
#     rate = 1.0
#     eqstrain_local = eqstrain_local
#     length_scale = ${fparse l}
#     kappa_i = ${fparse kappa_i}
#     c0 = ${fparse c0}
#   []
#   [diffusion_nonlocal]
#     type = CoefDiffusion
#     variable = nonlocal_eqstrain
#     coef = ${fparse 1.0}
#   []
#   [reaction_local]
#     type = CoupledElkLocalEqstrainForce
#     variable = nonlocal_eqstrain
#     eqstrain_local = eqstrain_local
#     length_scale = ${fparse l}
#     kappa_i = ${fparse kappa_i}
#     c0 = ${fparse c0}
#   []    
# []

# [Preconditioning]
#     [smp]
#       type = SMP
#       full = true
#     []
# []

[Executioner]
  type = Transient

  solve_type = JFNK

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  # petsc_options_value = '101                asm      lu'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  # petsc_options_value = ' lu       mumps       100'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  nl_max_its = 50
[]

[Outputs]
  exodus = true
  time_step_interval = 1000
  print_linear_residuals = false
  csv = true
[]