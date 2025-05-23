# Mandel's problem of consolodation of a drained medium
# Using the FullySaturatedDarcyBase and FullySaturatedMassTimeDerivative kernels
#
# A sample is in plane strain.
# -a <= x <= a
# -b <= y <= b
# It is squashed with constant force by impermeable, frictionless plattens on its top and bottom surfaces (at y=+/-b)
# Fluid is allowed to leak out from its sides (at x=+/-a)
# The porepressure within the sample is monitored.
#
# As is common in the literature, this is simulated by
# considering the quarter-sample, 0<=x<=a and 0<=y<=b, with
# impermeable, roller BCs at x=0 and y=0 and y=b.
# Porepressure is fixed at zero on x=a.
# Porepressure and displacement are initialised to zero.
# Then the top (y=b) is moved downwards with prescribed velocity,
# so that the total force that is inducing this downwards velocity
# is fixed.  The velocity is worked out by solving Mandel's problem
# analytically, and the total force is monitored in the simulation
# to check that it indeed remains constant.
#
# Here are the problem's parameters, and their values:
# Soil width.  a = 1
# Soil height.  b = 0.1
# Soil's Lame lambda.  la = 0.5
# Soil's Lame mu, which is also the Soil's shear modulus.  mu = G = 0.75
# Soil bulk modulus.  K = la + 2*mu/3 = 1
# Drained Poisson ratio.  nu = (3K - 2G)/(6K + 2G) = 0.2
# Soil bulk compliance.  1/K = 1
# Fluid bulk modulus.  Kf = 8
# Fluid bulk compliance.  1/Kf = 0.125
# Soil initial porosity.  phi0 = 0.1
# Biot coefficient.  alpha = 0.6
# Biot modulus.  M = 1/(phi0/Kf + (alpha - phi0)(1 - alpha)/K) = 4.705882
# Undrained bulk modulus. Ku = K + alpha^2*M = 2.694118
# Undrained Poisson ratio.  nuu = (3Ku - 2G)/(6Ku + 2G) = 0.372627
# Skempton coefficient.  B = alpha*M/Ku = 1.048035
# Fluid mobility (soil permeability/fluid viscosity).  k = 1.5
# Consolidation coefficient.  c = 2*k*B^2*G*(1-nu)*(1+nuu)^2/9/(1-nuu)/(nuu-nu) = 3.821656
# Normal stress on top.  F = 1
#
# The solution for porepressure and displacements is given in
# AHD Cheng and E Detournay "A direct boundary element method for plane strain poroelasticity" International Journal of Numerical and Analytical Methods in Geomechanics 12 (1988) 551-572.
# The solution involves complicated infinite series, so I shall not write it here

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator
  block = 0
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
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
    [roller_yminuu]
    type = DirichletBC
    variable = disp_y
    value = 1
    boundary = 'top'
  []
  [xmax_drained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = right
  []
  [top_velocity]
    type = NeumannBC
    variable = porepressure
    value = 0
    boundary = left
  []
[]

[Kernels]
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    variable = disp_y
    component = 1
  []
  [mass0]
    type = PorousFlowFullySaturatedMassTimeDerivative
    biot_coefficient = 0.6
    coupling_type = HydroMechanical
    variable = porepressure
  []
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = porepressure
    gravity = '0 0 0'
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 8
    density0 = 1
    thermal_expansion = 0
    viscosity = 1
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
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
  [eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
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
  end_time = 0.0001
  dt = 0.0001
[]

[Outputs]
    exodus = true
[]