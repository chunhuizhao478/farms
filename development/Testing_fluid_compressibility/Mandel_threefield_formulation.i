

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
  elem_type = QUAD9
[]

[GlobalParams]
    displacements = 'disp_x disp_y' 
   # fluid_vel = 'fluid_vel_x fluid_vel_y'
   # porepressure = 'p'
[]

[Variables]
  [disp_x]
        order = SECOND
        family = LAGRANGE
  []
  [disp_y]
        order = SECOND
        family = LAGRANGE
  []
  [fluid_vel_x]
        order = SECOND
        family = LAGRANGE
  []
  [fluid_vel_y]
        order = SECOND
        family = LAGRANGE
  []
  [p]
        order = FIRST
        family = LAGRANGE
  []
[]

[BCs]
  [roller_xmin]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'
  []
  [roller_ymin]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'
  []
  [xmax_drained]
    type = DirichletBC
    variable = fluid_vel_x
    value = 0
    boundary = left
  []
  [xmax_drained_2]
    type = DirichletBC
    variable = fluid_vel_y
    value = 0
    boundary = top
  []
  [xmax_drained_3]
    type = DirichletBC
    variable = fluid_vel_y
    value = 0
    boundary = bottom
  []
  [top_velocity]
    type = FunctionDirichletBC
    variable = disp_y
    function = top_velocity
    boundary = top
  []
[]

[Functions]
  [top_velocity]
    type = PiecewiseLinear
    x = '0 0.002 0.006   0.014   0.03    0.046   0.062   0.078   0.094   0.11    0.126   0.142   0.158   0.174   0.19 0.206 0.222 0.238 0.254 0.27 0.286 0.302 0.318 0.334 0.35 0.366 0.382 0.398 0.414 0.43 0.446 0.462 0.478 0.494 0.51 0.526 0.542 0.558 0.574 0.59 0.606 0.622 0.638 0.654 0.67 0.686 0.702'
    y = '-0.041824842    -0.042730269    -0.043412712    -0.04428867     -0.045509181    -0.04645965     -0.047268246 -0.047974749      -0.048597109     -0.0491467  -0.049632388     -0.050061697      -0.050441198     -0.050776675     -0.051073238      -0.0513354 -0.051567152      -0.051772022     -0.051953128 -0.052113227 -0.052254754 -0.052379865 -0.052490464 -0.052588233 -0.052674662 -0.052751065 -0.052818606 -0.052878312 -0.052931093 -0.052977751 -0.053018997 -0.053055459 -0.053087691 -0.053116185 -0.053141373 -0.05316364 -0.053183324 -0.053200724 -0.053216106 -0.053229704 -0.053241725 -0.053252351 -0.053261745 -0.053270049 -0.053277389 -0.053283879 -0.053289615'
  []
[]

[Kernels]
   [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false     
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_x
        porepressure = p
        component = 0
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
    [../]
    [./darcyflow_x]
        type = DarcyFlow
        variable = fluid_vel_x
    [../]
    [./darcyflow_y]
        type = DarcyFlow
        variable = fluid_vel_y
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_x
        porepressure = p
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_y
        porepressure = p
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p
        displacements = 'disp_x disp_y'
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p
        u = fluid_vel_x
        v = fluid_vel_y
        pressure = p
    [../]
[]


[Materials]
    [elasticity_tensor]
      type = ComputeElasticityTensor
      C_ijkl = '0.5 0.75'
      # bulk modulus is lambda + 2*mu/3 = 0.5 + 2*0.75/3 = 1
      fill_method = symmetric_isotropic
    []
    [strain]
      type = ComputeSmallStrain
    []
    [stress]
      type = ComputeLinearElasticStress
    []
    [./porosity]
        type = GenericConstantMaterial
        prop_names = porosity
        prop_values = 0.1
    [../]
    [./hydconductivity]
        type = GenericConstantMaterial
        prop_names = hydconductivity
        prop_values = 1.5
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.6
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 4.7058823529
    [../]
    [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
    [../]
[]

[Preconditioning]
    [./smp]
        type = SMP
        full = true
       # petsc_options_iname = '-snes_test_jacobian  -pc_factor_mat_solver_type -pc_factor_shift_type'
        #petsc_options_value = 'lu umfpack NONZERO'
    [../]
[]

[Executioner]
    type = Transient
    solve_type = Newton
    dt = 0.001
    start_time = 0
    end_time = 0.7
   # automatic_scaling = true  
[]

[Outputs]
    exodus = true
[]