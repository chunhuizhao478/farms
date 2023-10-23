[Mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 80
    xmin = 0
    xmax = 0.1
    ymin = 0
    ymax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator
  block = 0
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [porepressure]
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
  [plane_strain]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'back front'
  []
  [xmax_drained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = right
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

[AuxVariables]
  [tot_force]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [tot_force]
    type = ParsedAux
    coupled_variables = 'stress_yy porepressure'
    execute_on = timestep_end
    variable = tot_force
    expression = '-stress_yy+0.6*porepressure'
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0.0
    bulk_modulus = 8.0
    viscosity = 1.0
    density0 = 1.0
  []
[]

[PorousFlowBasicTHM]
  coupling_type = HydroMechanical
  displacements = 'disp_x disp_y disp_z'
  multiply_by_density = false
  porepressure = porepressure
  biot_coefficient = 0.6
  gravity = '0 0 0'
  fp = the_simple_fluid
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
  [porosity]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.1
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.6
    solid_bulk_compliance = 1
    fluid_bulk_modulus = 8
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.5 0 0   0 1.5 0   0 0 1.5'
  []
[]

[Postprocessors]
  [p0]
    type = PointValue
    outputs = csv
    point = '0.0 0 0'
    variable = porepressure
  []
  [p1]
    type = PointValue
    outputs = csv
    point = '0.1 0 0'
    variable = porepressure
  []
  [p2]
    type = PointValue
    outputs = csv
    point = '0.2 0 0'
    variable = porepressure
  []
  [p3]
    type = PointValue
    outputs = csv
    point = '0.3 0 0'
    variable = porepressure
  []
  [p4]
    type = PointValue
    outputs = csv
    point = '0.4 0 0'
    variable = porepressure
  []
  [p5]
    type = PointValue
    outputs = csv
    point = '0.5 0 0'
    variable = porepressure
  []
  [p6]
    type = PointValue
    outputs = csv
    point = '0.6 0 0'
    variable = porepressure
  []
  [p7]
    type = PointValue
    outputs = csv
    point = '0.7 0 0'
    variable = porepressure
  []
  [p8]
    type = PointValue
    outputs = csv
    point = '0.8 0 0'
    variable = porepressure
  []
  [p9]
    type = PointValue
    outputs = csv
    point = '0.9 0 0'
    variable = porepressure
  []
  [p99]
    type = PointValue
    outputs = csv
    point = '1 0 0'
    variable = porepressure
  []
  [xdisp]
    type = PointValue
    outputs = csv
    point = '1 0.1 0'
    variable = disp_x
  []
  [ydisp]
    type = PointValue
    outputs = csv
    point = '1 0.1 0'
    variable = disp_y
  []
  [total_downwards_force]
     type = ElementAverageValue
     outputs = csv
     variable = tot_force
  []
  [dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = if(0.15*t<0.01,0.15*t,0.01)
  []
[]

[Preconditioning]
  [andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres asm lu 1E-14 1E-10 10000'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = 0
  end_time = 0.7
  [TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 0.001
  []
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = mandel_basicthm
  [csv]
    interval = 3
    type = CSV
  []
[]