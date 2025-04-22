[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole.msh'
  []
  [markall_sideset]
    type = SideSetsAroundSubdomainGenerator
    input = msh
    block = 3
    new_boundary = 'allsidesets'
  []
[]

  [Variables]
    [nonlocal_eqstrain_sub]
        order = FIRST
        family = MONOMIAL
    []
  []

  [AuxVariables]
    [crack_damage_aux_sub]
      order = FIRST
      family = MONOMIAL
    []
    [eqstrain_local_aux_sub]
      order = FIRST
      family = MONOMIAL
    []
    [elastic_energy_aux_sub]
      order = CONSTANT
      family = MONOMIAL
    []
  []

  [Kernels]
    # gradient based nonlocal averaging
    [react_nonlocal]
        type = Reaction
        variable = nonlocal_eqstrain_sub
        rate = 1.0
    []
    [diffusion_nonlocal]
        type = CoefDiffusion
        variable = nonlocal_eqstrain_sub
        coef = 1e-9 #1e-4
    []
    [reaction_local]
        type = CoupledElkLocalEqstrainForce
        eqstrain_local = eqstrain_local_aux_sub
        variable = nonlocal_eqstrain_sub
    []    
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
  []

  [BCs]
    [natural_bc]
      type = NeumannBC
      variable = nonlocal_eqstrain_sub
      boundary = allsidesets
      value = 0      
    []
  []

  [Adaptivity]
    max_h_level = 2
    marker = 'combo'
    [Indicators]
        [error]
          type = GradientJumpIndicator
            variable = crack_damage_aux_sub
        []
    []
    [Markers]
        [./combo]
            type = ComboMarker
            markers = 'error_marker elastic_energy_marker'
        [../]
        [./error_marker]
            type = ErrorFractionMarker
            indicator = error
            refine = 0.9
        [../]
        # peak energy: 144
        [./elastic_energy_marker]
            type = ValueThresholdMarker
            variable = elastic_energy_aux_sub
            refine = 100
        []         
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = NEWTON #in smeared cracking w/o full jacobian, this is much more efficient than NEWTON
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol' 
    petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9'
    # automatic_scaling = true
    # line_search = 'bt'
    # num_steps = 1
    l_max_its = 50
    nl_max_its = 100
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    # start_time = 0.0
    # end_time = 100
    # dt = 0.0001
    # verbose = true
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []