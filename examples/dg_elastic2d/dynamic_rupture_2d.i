# Minimal 2D DG elastic dynamic rupture test (plane strain) using new generalized kernels.
# Mesh: square domain with fault (y=0) split into two blocks to allow interface DG fault kernel.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
    xmin = -10000
    xmax =  10000
    ymin = -5000
    ymax =  5000
  []
  [upper]
    type = ParsedSubdomainMeshGenerator
    input = gen
    combinatorial_geometry = 'y>0'
    block_id = 1
  []
  [lower]
    type = ParsedSubdomainMeshGenerator
    input = upper
    combinatorial_geometry = 'y<0'
    block_id = 2
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = lower
    split_interface = true
    block_pairs = '1 2'
  []
[]

[GlobalParams]
  lambda_name = lambda
  mu_name = mu
  rho_name = rho
[]

[Variables]
  [./ux]
    order = FIRST
    family = LAGRANGE
  [../]
  [./uy]
    order = FIRST
    family = LAGRANGE
  [../]
  # Stresses (plane strain subset) sxx, syy, szz, sxy retained; szz evolves for plane strain closure
  [./sxx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./syy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./szz]
    order = FIRST
    family = LAGRANGE
  [../]
  [./sxy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./slip]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Materials]
  [./elastic]
    type = GenericConstantMaterial
    prop_names = 'lambda mu rho'
    prop_values = '32.04e9 32.04e9 2670'
  [../]
[]

# Time derivative (mass) terms for velocities with density
[Kernels]
  # Velocity equations: rho * du/dt handled by RhoTimeDerivative + flux from stresses
  [./rho_dt_ux]
    type = RhoTimeDerivative
    variable = ux
    rho_name = rho
  [../]
  [./rho_dt_uy]
    type = RhoTimeDerivative
    variable = uy
    rho_name = rho
  [../]

  # Velocity flux kernels (divergence of stress columns)
  [./vel_flux_ux]
    type = ElasticVelocityFlux3D
    variable = ux
    component = ux
    sxx = sxx
    syy = syy
    szz = szz
    sxy = sxy
    sxz = sxy   # dummy reuse (2D) not used in residual
    syz = sxy   # dummy reuse (2D)
  [../]
  [./vel_flux_uy]
    type = ElasticVelocityFlux3D
    variable = uy
    component = uy
    sxx = sxx
    syy = syy
    szz = szz
    sxy = sxy
    sxz = sxy
    syz = sxy
  [../]

  # Stress evolution kernels from velocities
  [./flux_sxx]
    type = ElasticStressFlux3D
    variable = sxx
    component = sxx
    ux = ux
    uy = uy
    uz = uy # 2D placeholder
    lambda_name = lambda
    mu_name = mu
  [../]
  [./flux_syy]
    type = ElasticStressFlux3D
    variable = syy
    component = syy
    ux = ux
    uy = uy
    uz = uy
    lambda_name = lambda
    mu_name = mu
  [../]
  [./flux_szz]
    type = ElasticStressFlux3D
    variable = szz
    component = szz
    ux = ux
    uy = uy
    uz = uy
    lambda_name = lambda
    mu_name = mu
  [../]
  [./flux_sxy]
    type = ElasticStressFlux3D
    variable = sxy
    component = sxy
    ux = ux
    uy = uy
    uz = uy
    lambda_name = lambda
    mu_name = mu
  [../]
[]

# Fault interface DG kernel with slip-weakening (placeholder usage). Provide minus/plus traces of fields.
[DGKernels]
  [./fault]
    type = FaultSlipWeakeningDG3D
    variable = ux   # one DGKernel per primal variable is typical; replicate for others if needed
    ux = ux
    uy = uy
    uz = uy          # 2D placeholder
    sxx = sxx
    syy = syy
    szz = szz
    sxy = sxy
    sxz = sxy
    syz = sxy
    slip = slip
    lambda_name = lambda
    mu_name = mu
    rho_name = rho
    mu_s = 0.677
    mu_d = 0.525
    Dc = 0.4
    boundary = '1_2'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.0005
  end_time = 0.05
  [TimeIntegrator]
    type = ExplicitEuler
  []
[]

[Outputs]
  exodus = true
  csv = true
  interval = 10
[]
