[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../mesh/cylinderlocal_refined.msh'
  []
  displacements = 'disp_x disp_y disp_z'
[]
  
  [GlobalParams]
    displacements = 'disp_x disp_y disp_z'
  []

  [Variables]
    [disp_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
    []
    [disp_z]
        order = FIRST
        family = LAGRANGE
    []
  []

  [AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [vel_z]
        order = FIRST
        family = LAGRANGE
    []
    [accel_z]
        order = FIRST
        family = LAGRANGE
    []
    [./strength]
      order = CONSTANT
      family = MONOMIAL
    [../]
  []
  
  [Functions]
    # [func_bc]
    #     type = PiecewiseConstant
    #     data_file = "../pulsepower/pressuredata_2.csv"
    #     format = columns
    # []
    [func_tri_pulse]
      type = ParsedFunction
      expression = '-10e12 * t + 80e6'
    []
  []

  [Kernels]
    [dispkernel_x]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
    []
    [dispkernel_y]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
    []
    [dispkernel_z]
        type = ADDynamicStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
        use_displaced_mesh = true
    []
    [inertia_x]
        type = ADInertialForce
        variable = disp_x
        velocity = vel_x
        acceleration = accel_x
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_y]
        type = ADInertialForce
        variable = disp_y
        velocity = vel_y
        acceleration = accel_y
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
    [inertia_z]
        type = ADInertialForce
        variable = disp_z
        velocity = vel_z
        acceleration = accel_z
        beta = 0.25
        gamma = 0.5
        use_displaced_mesh = true
    []
  []

  [AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = timestep_end
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = timestep_end
    []
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = timestep_end
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = timestep_end
    []
  []
  
  [BCs]
    [./Pressure]
        [pulsepower_load]
            boundary = 4
            function = func_tri_pulse
            displacements = 'disp_x disp_y'
        []
        # [static_pressure_outer]
        #   boundary = 3
        #   factor = 10e6
        #   displacements = 'disp_x disp_y'
        # []
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ADComputeIsotropicElasticityTensor
      youngs_modulus = 40e9
      poissons_ratio = 0.25
    [../]
    [./elastic_stress]
      type = ADComputeSmearedCrackingStress
      cracking_stress = strength
      softening_models = abrupt_softening
      cracked_elasticity_type = FULL
      max_cracks = 1
      output_properties = 'crack_damage stress'
      outputs = exodus
    [../]
    [strain]
        type = ADComputeFiniteStrain
        displacements = 'disp_x disp_y disp_z'
    []
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = '2600'
    []
    [./abrupt_softening]
      type = ADAbruptSoftening
    [../]
    [./exponential_softening]
      type = ADExponentialSoftening
    [../]
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []

  [ICs]
    [./strength]
      type = VolumeWeightedWeibull
      variable = strength
      reference_volume = 1e-7
      weibull_modulus = 12.0
      median = 30e6
    [../]
  []
  
  [Executioner]
    type = Transient
    solve_type = Newton
    petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    petsc_options_value = '101                asm      lu'
  
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 10e-6
    dt = 1e-8
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
  []
  
  # [Outputs]
  #   exodus = true
  #   interval = 5
  # []