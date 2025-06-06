[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../../meshfile/fieldscale_test2.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
[]

[Adaptivity]
  max_h_level = 5
  marker = 'combo'
  cycles_per_step = 1
  [Markers]
      [./combo]
        type = FarmsComboMarker
        markers = 'damage_marker strain_energy_marker'
        meshsize_marker = 'meshsize_marker'
        block = 2
      [../]
      [damage_marker]
        type = ValueThresholdMarker
        variable = d
        refine = 0.5
        block = 2
      []
      [strain_energy_marker]
        type = ValueThresholdMarker
        variable = psie_active
        refine = '${fparse 1.0*3/8*Gc_const/l}'
        block = 2
      []   
      # if mesh_size > dxmin, refine
      # if mesh_size < dxmin/100, coarsen (which never happens)
      # otherwise, do nothing
      [meshsize_marker]
        type = ValueThresholdMarker
        variable = mesh_size
        refine = '${fparse 2 * dx_min}'
        coarsen = '${fparse dx_min/100}'
        third_state = DO_NOTHING
        block = 2
      [] 
  []
[]

[Variables]
  [d]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [Gc_var]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []
[]

[Kernels]
  [diff]
    type = ADPFFDiffusion #
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l'
    prop_values = '${l}'
  []
  [Gc_var]
    type = ADParsedMaterial
    property_name = Gc
    coupled_variables = 'Gc_var'
    expression = 'Gc_var'
    # outputs = exodus
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [crack_geometric] #AT1 model
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
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
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  # petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
  petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'

  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  # time_step_interval = 40
  print_linear_residuals = false
[]

[Distributions]
  #typically for granite
  #Shape Parameter (k): 5 to 15, commonly around 8 to 12.
  #Scale Parameter (λ): 5 to 30 MPa, commonly around 10 to 20 MPa.
  [weibull]
    type = Weibull
    shape = 12.0 #k
    scale = ${Gc_const} #lambda
    location = 0 
  []
[] 

[ICs]
  [./gc_var]
    type =  RandomIC
    variable = Gc_var
    distribution = weibull
  []
[]