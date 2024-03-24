[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 721
        ny = 721
        xmin = -360
        xmax = 360
        ymin = -720
        ymax = 0
    []
[]

[GlobalParams]
    q = 0.2
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
[]

[AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
[]

[Physics/SolidMechanics/QuasiStatic]
    displacements = 'disp_x disp_y'
    [./main]
      strain = SMALL
      add_variables = true
      displacements = 'disp_x disp_y'
    [../]
[]

[Kernels]
    [./inertia_x]
        type = InertialForce
        use_displaced_mesh = true
        variable = disp_x
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = true
        variable = disp_y
    []
    [./damping_x]
        type = StiffPropDamping
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
    []
    [./damping_y]
        type = StiffPropDamping
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1       
    []
[]

[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[Materials]
    [stress_medium]
        type = ComputeLinearElasticStress
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        shear_modulus = 32.038e9
        poissons_ratio = 0.25
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
[]

[Functions]
    [./sbi_inputdata_vx]
        type = PiecewiseMultilinear
        data_file = ../events_data_piecewisemultilinear/events_data_vx_num13.txt
    [../]
    [./sbi_inputdata_vy]
        type = PiecewiseMultilinear
        data_file = ../events_data_piecewisemultilinear/events_data_vy_num13.txt
    [../]
[]

[Times]
    [file]
        type = CSVFileTimes
        files = ../events_data_piecewisemultilinear/events_data_time_num13.csv
    []
[]
  
[Executioner]
    type = Transient
    dt = 6.94e-5
    num_steps = 4000
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []    
[]

[Outputs]
    exodus = true
    interval = 10
    show = 'vel_x vel_y'
    [checkpoints]
      type = Checkpoint
      num_files = 1
      execute_on = 'FINAL'
    []
[]

[BCs]
    #assume we are doing sbi_bottom, the top boundary is sbi data, other boundaries are absorbing boundary conditions
    #sbi input data top boundaries
    [./sbi_inputdata_top_x]
        type = PresetVelocity
        variable = disp_x
        boundary = top
        function = sbi_inputdata_vx
    []
    [./sbi_inputdata_top_y]
        type = PresetVelocity
        variable = disp_y
        boundary = top
        function = sbi_inputdata_vy
    []
    #absorbing bcs for bottom, left, right
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 5999.81
        shear_wave_speed = 3464
        boundary = right
    []
[]