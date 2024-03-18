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
        data_file = ../events_data_piecewisemultilinear/events_data_vx_num1.txt
    [../]
    [./sbi_inputdata_vy]
        type = PiecewiseMultilinear
        data_file = ../events_data_piecewisemultilinear/events_data_vy_num1.txt
    [../]
[]

[Times]
    [file]
        type = CSVFileTimes
        files = ../events_data_piecewisemultilinear/events_data_time_num1.csv
    []
[]
  
[Executioner]
    type = Transient

    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []    
  
    [TimeStepper]
      type = TimeSequenceFromTimes
      times = file
    []
[]

[Outputs]
    exodus = true
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
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6095
        shear_wave_speed = 3129
        boundary = right
    []
[]