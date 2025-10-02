# 2D DG elastic test with Godunov flux on all internal faces except the fault, where slip-weakening acts.

[Mesh]
  # Base structured grid covering full domain
  [gen]
    type = GeneratedMeshGenerator
    dim = 2                    # 2D plane-strain test
    nx = 40                    # elements in x (adjust for resolution study)
    ny = 40                    # elements in y
    xmin = -10000 
    xmax =  10000
    ymin = -5000  
    ymax =  5000
  []
  # Tag upper half-space as block 1 (y>0)
  [upper]
    type = ParsedSubdomainMeshGenerator
    input = gen
    combinatorial_geometry = 'y>0'
    block_id = 1
  []
  # Tag lower half-space as block 2 (y<0)
  [lower]
    type = ParsedSubdomainMeshGenerator
    input = upper
    combinatorial_geometry = 'y<0'
    block_id = 2
  []
  # Split interface between block 1 & 2 so DGKernels can see a distinct interior boundary (fault)
  [./iface]
    type=SideSetsBetweenSubdomainsGenerator
    input=lower
    primary_block='1'
    paired_block='2'
    new_boundary='Block1_Block2'
  [../]
[]

[GlobalParams]
  # Material property name aliases used by kernels / DG kernels
  lambda_name = lambda
  mu_name = mu
  rho_name = rho
[]

[Variables]
  # Primary velocity components
  [./ux] order=FIRST family=LAGRANGE [../]
  [./uy] order=FIRST family=LAGRANGE [../]
  # Stress tensor components (plane strain keeps szz to enforce constitutive constraint)
  [./sxx] order=FIRST family=LAGRANGE [../]
  [./syy] order=FIRST family=LAGRANGE [../]
  [./szz] order=FIRST family=LAGRANGE [../]
  [./sxy] order=FIRST family=LAGRANGE [../]
[]

[Materials]
  # Uniform elastic and density properties (can later make spatially variable).
  [./elastic]
    type = GenericConstantMaterial
    prop_names = 'lambda mu rho'
    prop_values = '32.04e9 32.04e9 2670'
  [../]
  # Fault interface material computing slip-weakening fluxes & diagnostics
  [./fault_mat]
    type = FaultSlipWeakeningDG3DMaterial
    boundary = Block1_Block2
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    mu_s=0.677 mu_d=0.525 Dc=0.4
    tau0_t1=0.0 tau0_t2=0.0 sigma0n=0.0
    mu_s_aux=mus tau0_t1_aux=tau0_t1_field sigma0n_aux=sigma0n_field
    boundary = 'Block1_Block2'
  [../]
[]

[Functions]
  [func_static_friction_coeff_mus]
    type = PiecewiseConstant
    axis=x
    x = '-1000e3 -15e3 15e3'
    y = '10000 0.677 10000.0'
    direction = left
  []
  [func_initial_strike_shear_stress]
    type = PiecewiseConstant
    axis=x
    x = '-1000e3 -9.0e3 -6.0e3 -1.5e3  1.5e3  6.0e3  9.0e3'
    y = ' 70.0e6 78.0e6 70.0e6 81.6e6 70.0e6 62.0e6 70.0e6'
  []
  [func_initial_normal_stress]
    type = ConstantFunction
    value = 120e6
  []
[]

[Kernels]
  # Density-weighted time derivative (ρ ∂t u) for ux
  [rho_dt_ux] type=RhoTimeDerivative variable=ux rho_name=rho [../]
  # Density-weighted time derivative (ρ ∂t u) for uy
  [rho_dt_uy] type=RhoTimeDerivative variable=uy rho_name=rho [../]
  # Stress time derivatives (unit coefficients) ensuring evolution equations, ∂t σ_ij
  [dt_sxx] type=TimeDerivative variable=sxx [../]
  [dt_syy] type=TimeDerivative variable=syy [../]
  [dt_szz] type=TimeDerivative variable=szz [../]
  [dt_sxy] type=TimeDerivative variable=sxy [../]
  # Velocity flux (divergence of stress column) for ux equation
  [vel_flux_ux]
    type = ElasticVelocityFlux3D
    variable = ux
    component = ux
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
  [../]
  # Velocity flux for uy equation
  [vel_flux_uy]
    type = ElasticVelocityFlux3D
    variable = uy
    component = uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
  [../]
  # Stress evolution fluxes from velocity field gradients
  [flux_sxx]
    type = ElasticStressFlux3D
    variable = sxx
    component = sxx
    ux=ux uy=uy uz=uy  # uz placeholder in 2D
  [../]
  [flux_syy]
    type = ElasticStressFlux3D
    variable = syy
    component = syy
    ux=ux uy=uy uz=uy
  [../]
  [flux_szz]
    type = ElasticStressFlux3D
    variable = szz
    component = szz
    ux=ux uy=uy uz=uy
  [../]
  [flux_sxy]
    type = ElasticStressFlux3D
    variable = sxy
    component = sxy
    ux=ux uy=uy uz=uy
  [../]
[]

[DGKernels]
  # Godunov numerical flux (exact elastic Riemann solver) applied to all internal faces
  # except those enumerated in exclude_boundaries (the fault).
  [godunov_ux]
    type = ElasticGodunovDGFlux3D
    variable = ux
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]
  [godunov_uy]
    type = ElasticGodunovDGFlux3D
    variable = uy
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]
  # Godunov fluxes for stress components (propagate stresses via velocity jumps)
  [godunov_sxx]
    type = ElasticGodunovDGFlux3D
    variable = sxx
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]
  [godunov_syy]
    type = ElasticGodunovDGFlux3D
    variable = syy
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]
  [godunov_szz]
    type = ElasticGodunovDGFlux3D
    variable = szz
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]
  [godunov_sxy]
    type = ElasticGodunovDGFlux3D
    variable = sxy
    ux=ux uy=uy uz=uy
    sxx=sxx syy=syy szz=szz sxy=sxy sxz=sxy syz=sxy
    lambda_name=lambda mu_name=mu rho_name=rho
    skip_block_pair = '1 2'
  [../]

  # Fault interface DG kernels applying slip-weakening friction law on boundary 1_2
  [fault_ux]
    type = FaultSlipWeakeningDG3D
    variable = ux
    boundary = 'Block1_Block2'
  # Fluxes now supplied purely via material properties
  [../]
  [fault_uy]
    type = FaultSlipWeakeningDG3D
    variable = uy
    boundary = 'Block1_Block2'
  # Fluxes supplied by fault_mat
  [../]
  [fault_sxx]
    type = FaultSlipWeakeningDG3D
    variable = sxx
    boundary = 'Block1_Block2'
  # Material provides flux
  [../]
  [fault_syy]
    type = FaultSlipWeakeningDG3D
    variable = syy
    boundary = 'Block1_Block2'
  # Material provides flux
  [../]
  [fault_szz]
    type = FaultSlipWeakeningDG3D
    variable = szz
    boundary = 'Block1_Block2'
  # Material provides flux
  [../]
  [fault_sxy]
    type = FaultSlipWeakeningDG3D
    variable = sxy
    boundary = 'Block1_Block2'
  # Material provides flux
  [../]
[]

## Removed FaultGodunovStateUO: superseded by fault interface material

[AuxVariables]
  [./mus]  order=FIRST family=LAGRANGE [../]
  [./tau0_t1_field] order=FIRST family=LAGRANGE [../]
  [./sigma0n_field] order=FIRST family=LAGRANGE [../]
[]

[AuxKernels]
  [mus_init]
    type=FunctionAux
    variable=mus
    function=func_static_friction_coeff_mus
    boundary='Block1_Block2'
  []
  [tau0t1_init]
    type=FunctionAux
    variable=tau0_t1_field
    function=func_initial_strike_shear_stress
    boundary='Block1_Block2'
  []
  [sigma0n_init]
    type=FunctionAux
    variable=sigma0n_field
    function=func_initial_normal_stress
    boundary='Block1_Block2'
  []
[]

[Executioner]
  type = Transient
  dt = 0.01
  end_time = 12
  # num_steps = 1
  [TimeIntegrator]
    type = ExplicitEuler
  []
[]


