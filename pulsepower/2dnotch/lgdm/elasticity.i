E = 50e9
nu = 0.373
ft = 137e6 ##computed from pf
# Gc_const = 100
density = 2600
# dx_min = 5e-5

h_modulus = '${fparse 1e-9 * E}'

K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  1e-4 
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/density)}'

#gradient activity parameters
kappa_i = ${fparse ft / E}
c0 = 1e-12 #minimum value of the gradient activity parameter for the equivalent strain

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0
#----------------------------------------------------#

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = nonlocal_subapp.i
    cli_args = 'l=${l};kappa_i=${kappa_i};c0=${c0}'
    execute_on = 'TIMESTEP_BEGIN'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = nonlocal_eqstrain
    source_variable = nonlocal_eqstrain
    execute_on = 'TIMESTEP_BEGIN'
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'eqstrain_local crack_damage_aux'
    source_variable = 'eqstrain_local crack_damage_aux'
    execute_on = 'TIMESTEP_BEGIN'
  []
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
  [./strength]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${fparse ft}
  [../]
  [crack_damage_aux]
    order = FIRST
    family = MONOMIAL
  []
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    order = CONSTANT
    family = MONOMIAL
  []
  [crack_damage_initial]
    family = LAGRANGE
    order = FIRST
  []
  [nonlocal_eqstrain]
    order = FIRST
    family = LAGRANGE
  [] 
  [eqstrain_local]
    family = MONOMIAL
    order = CONSTANT
  []
  [accel_x]
  []
  [accel_y]
  []
  [vel_x]
  []
  [vel_y]
  []
  [vel_z]
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
  #
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  #
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  #get pulse load aux
  # [get_pulse_load_aux]
  #   type = FunctionAux 
  #   variable = pulse_load_aux
  #   function = func_tri_pulse
  #   execute_on = timestep_end
  # []
  #mesh size aux
  [./max]
    type = ElementLengthAux
    variable = mesh_size
    method = max
    execute_on = TIMESTEP_BEGIN
  [../]
  #damage
  # [define_initial_damage_block1]
  #   type = ConstantAux
  #   variable = crack_damage_initial
  #   value = 0
  #   block = 1
  #   execute_on = INITIAL
  # []
  # [define_initial_damage_block0]
  #   type = ConstantAux
  #   variable = crack_damage_initial
  #   value = 0
  #   block = '4 5'
  #   execute_on = INITIAL
  # []
  #get eqstrain_local
  [eqstrain_local_aux]
    type = MaterialRealAux
    variable = eqstrain_local
    property = eqstrain_local
    execute_on = 'INITIAL NONLINEAR TIMESTEP_END'
  []
  #get crack damage aux
  [crack_damage_aux]
    type = MaterialRealAux
    variable = crack_damage_aux
    property = crack_damage
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
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = fy
  []  
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
  [./strain]
    type = ComputeSmallStrain
  []
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [./elastic_stress]
    type = FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain
    nonlocal_eqstrain = nonlocal_eqstrain
    paramA = 0.99
    paramB = 500
    cracking_stress = strength
    initial_crack_damage = crack_damage_initial
    output_properties = 'elastic_strain psie_active strain_increment'
    h = ${h_modulus}
    outputs = exodus
  [../]
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${density}
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

  solve_type = 'NEWTON'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  # petsc_options_value = '101                asm      lu'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  # petsc_options_value = ' lu       mumps       100'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 30

  dt = 1e-8
  end_time = 1e10

  # fixed_point_max_its = 10
  # accept_on_max_fixed_point_iteration = false
  # fixed_point_rel_tol = 1e-6
  # fixed_point_abs_tol = 1e-8

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
  [./TimeIntegrator]
    type = ImplicitEuler
  [../]
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
  time_step_interval = 50
  print_linear_residuals = false
  [csv]
    type = CSV
    execute_on = 'initial timestep_end'
    time_step_interval = 1
  []
[]