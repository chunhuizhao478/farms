E = 50e9
nu = 0.373
# ft = 25.5e6
Gc_const = 100
density = 2600
# dx_min = 5e-5

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  1e-5
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
# Cs = '${fparse sqrt(G/density)}'
# Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0 #match energy budget
#----------------------------------------------------#

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc_const=${Gc_const};l=${l}'
    execute_on = 'INITIAL TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = 'd dissipated_energy_density'
    source_variable = 'd dissipated_energy_density'
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active mesh_size'
  []
  # [pp_transfer_dissipated_energy_total]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = 'fracture'
  #   from_postprocessor = 'dissipated_energy_dynamic'
  #   to_postprocessor = 'dissipated_energy_dynamic'
  #   reduction_type = 'sum' #this should not have effect in a single app
  # []
  # [pp_transfer_dissipated_energy_first_step]
  #   type = MultiAppPostprocessorTransfer
  #   from_multi_app = 'fracture'
  #   from_postprocessor = 'dissipated_energy_first_step'
  #   to_postprocessor = 'dissipated_energy_first_step'
  #   reduction_type = 'sum' #this should not have effect in a single app
  # []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../mesh/2dphysicalnotch.msh'
    []
    displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = FIRST
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  #err measurement of active strain energy
  [eng_err]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  #
  [vel_z]
    family = LAGRANGE
    order = FIRST
  []
  #
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  #
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  #
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  #
  [dissipated_energy_density]
    family = MONOMIAL
    order = CONSTANT
  []
  #reaction force
  [fx]
  []
  [fy]
  []
  [fz]
  []
  [fdampx]
  []
  [fdampy]
  []
  [fdampz]
  []
[]

[AuxKernels]
  [./error_measure]
    type = ErrorPsiMeasure
    variable = eng_err
  [../]
  #
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  #
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # #get pulse load aux
  # [get_pulse_load_aux]
  #   type = FunctionAux 
  #   variable = pulse_load_aux
  #   function = func_tri_pulse
  #   execute_on = timestep_end
  # []
  #mesh size aux
  [mesh_size_aux]
    type = MeshSize
    variable = mesh_size
    execute_on = 'TIMESTEP_END'
  []
[]

[Functions]
  [func_dyn_load]
      type = ParsedFunction
      expression = '1 * t'
  []
[]

[Kernels]
  [solid_x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y'
  []
  [solid_y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    displacements = 'disp_x disp_y'
    save_in = fy
  []
  # [inertia_x]
  #   type = ADInertialForce
  #   variable = disp_x
  #   use_displaced_mesh = false
  #   beta = ${newmark_beta}
  #   gamma = ${newmark_gamma}
  #   velocity = vel_x
  #   acceleration = accel_x
  #   density = ${density}
  #   alpha = ${hht_alpha}
  # []
  # [inertia_y]
  #   type = ADInertialForce
  #   variable = disp_y
  #   use_displaced_mesh = false
  #   beta = ${newmark_beta}
  #   gamma = ${newmark_gamma}
  #   velocity = vel_y
  #   acceleration = accel_y
  #   density = ${density}
  #   alpha = ${hht_alpha}
  # []
[]

[BCs]
  [fix_bottom_x]
      type = DirichletBC
      variable = disp_x
      boundary = 2
      value = 0.0
  []
  [fix_bottom_y]
      type = DirichletBC
      variable = disp_y
      boundary = 2
      value = 0.0
  []
  [./load_top]
      type = FunctionDirichletBC
      variable = disp_y
      boundary = 1
      function = func_dyn_load
  [../]
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G density'
    prop_values = '${K} ${G} ${density}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active psie'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [strain]
    type = ADComputeSmallStrain
  []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '

  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  # petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true
  line_search = 'basic'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  # Add more iterations before failure
  nl_max_its = 30

  # dt = 0.5e-7
  end_time = 1e10

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  dt = 1e-8
  # [TimeStepper]
  #   type = FarmsIterationAdaptiveDT
  #   dt = 1e-8
  #   iteration_window = 0 #the adaptive time stepping happens at number of iterations <-> 'optimal_iterations plus/minus iteration_window'
  #   cutback_factor_at_failure = 0.5
  #   optimal_iterations = 20
  #   growth_factor = 1.25
  #   max_time_step_bound = 1e-7
  # []
  # [./TimeIntegrator]
  #   type = NewmarkBeta
  #   beta = ${newmark_beta}
  #   gamma = ${newmark_gamma}
  # [../]
[]

[Postprocessors]
  [Fx]
    type = NodalSum
    variable = fy
    boundary = 1
  []
[]

[Outputs]
  exodus = true
  time_step_interval = 1
  print_linear_residuals = false
  [csv]
    type = CSV
    execute_on = 'initial timestep_end'
    time_step_interval = 1
  []
[]