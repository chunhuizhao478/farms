[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../mesh/2dphysicalnotch.msh'
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
    R = 0.005
    eta = 2.5 #5
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

[Executioner]
  type = Transient

  solve_type = NEWTON

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