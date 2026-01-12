[Mesh]
  type = FileMesh
  file = bp1qd_fault.e
[]

[GlobalParams]
  mu      = 3.204e10       # Pa
  rho     = 2670.0         # kg/m^3
  cs      = 3464.0         # m/s
  eta     = ${mu}/(2*${cs})
  sigma_n = 50e6           # Pa
  Dc      = 8e-3           # m
  a0      = 0.010
  amax    = 0.025
  b0      = 0.015
  f0      = 0.6
  V0      = 1e-6           # m/s
  Vp      = 1e-9           # m/s
  H       = 15000.0        # m
  h       = 3000.0         # m
  Wf      = 40000.0        # m
[]

[Variables]
  [./u]       # antiplane displacement (out-of-plane)
    family = LAGRANGE
    order  = FIRST
  [../]
  [./theta]   # state var at fault nodes
    family = LAGRANGE
    order  = FIRST
  [../]
[]

[Kernels]
  # Antiplane elasticity: ∇·(mu ∇u) = 0
  [./diff]
    type     = Diffusion
    variable = u
  [../]
[]

[DGKernels]
  # Rate–state friction + radiation damping on internal fault
  [./rsf]
    type               = RSFrictionDG
    variable           = u
    boundary           = 'fault_left fault_right'   # apply on both faces
    sigma_n            = ${sigma_n}
    eta                = ${eta}
    Dc                 = ${Dc}
    f0                 = ${f0}
    V0                 = ${V0}
    a0                 = ${a0}
    amax               = ${amax}
    b0                 = ${b0}
    H                  = ${H}
    h                  = ${h}
    Wf                 = ${Wf}
    theta              = theta
    mu                 = ${mu}
  [../]
[]

[NodalKernels]
  # dθ/dt = 1 - V θ / Dc at fault nodes
  [./theta_aging_left]
    type        = ThetaAgingNodal
    variable    = theta
    boundary    = fault_left
    Dc          = ${Dc}
    u_primary   = u      # left node value
    u_neighbor  = u      # right node value
  [../]
  [./theta_aging_right]
    type        = ThetaAgingNodal
    variable    = theta
    boundary    = fault_right
    Dc          = ${Dc}
    u_primary   = u
    u_neighbor  = u
  [../]
[]

[ICs]
  [./u0]
    type     = ConstantIC
    variable = u
    value    = 0.0
  [../]
  # Aging-law steady state for V = Vp → theta = Dc/Vp
  [./theta0_left]
    type      = ConstantBC   # nodal init on boundary nodes works via BC at t=0
    variable  = theta
    boundary  = fault_left
    value     = ${Dc}/${Vp}
  [../]
  [./theta0_right]
    type      = ConstantBC
    variable  = theta
    boundary  = fault_right
    value     = ${Dc}/${Vp}
  [../]
[]

[Functions]
  [./left_bottom_disp]
    type  = ParsedFunction
    value = -0.5*Vp*t
  [../]
  [./right_bottom_disp]
    type  = ParsedFunction
    value =  0.5*Vp*t
  [../]
[]

[BCs]
  # Pin far sides to avoid rigid shift in u
  [./left_fix]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0.0
  [../]
  [./right_fix]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 0.0
  [../]

  # Drive bottom with opposite plate motions so relative V = Vp at depth
  [./bottom_left_drive]
    type      = FunctionDirichletBC
    variable  = u
    boundary  = bottom&block_1
    function  = left_bottom_disp
  [../]
  [./bottom_right_drive]
    type      = FunctionDirichletBC
    variable  = u
    boundary  = bottom&block_2
    function  = right_bottom_disp
  [../]
[]

[Executioner]
  type       = Transient
  scheme     = bdf2
  start_time = 0.0
  end_time   = 9.46e10          # ~3000 years
  dt         = 1.0
  dt_min     = 1e-6
  dt_max     = 3.15e7           # ~1 year
[]

[Postprocessors]
  # Example probes: you can add multiple depths
  [./u_r_7p5]
    type     = PointValue
    variable = u
    point    = '0.1 7500.0'
  [../]
  [./u_l_7p5]
    type     = PointValue
    variable = u
    point    = '-0.1 7500.0'
  [../]
[]

[Outputs]
  exodus = true
  [./csv]
    type       = CSV
    execute_on = 'timestep_end'
    outputs    = 'u_r_7p5 u_l_7p5'
  [../]
[]