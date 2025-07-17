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
  dim = 3
  nx = 40
  ny = 4
  nz = 4
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
  zmin = 0
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = dictator
  block = 0
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
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
    x = '0 0.002 0.006 0.014 0.03 0.046 0.062 0.078 0.094 0.11 0.126 0.142 0.158 0.174 0.19 0.206 0.222 0.238 0.254 0.27 0.286 0.302 0.318 0.334 0.35 0.366 0.382 0.398 0.414 0.43 0.446 0.462 0.478 0.494 0.51 0.526 0.542 0.558 0.574 0.59 0.606 0.622 0.638 0.654 0.67 0.686 0.702'
    y = '-0.0002091242 -0.0002136514 -0.0002170636 -0.0002214434 -0.0002275459 -0.0002322983 -0.0002363412 -0.0002398737 -0.0002429855 -0.0002457335 -0.0002481619 -0.0002503085 -0.0002522060 -0.0002538834 -0.0002553662 -0.0002566770 -0.0002578358 -0.0002588601 -0.0002597656 -0.0002605661 -0.0002612738 -0.0002618993 -0.0002624523 -0.0002629412 -0.0002633733 -0.0002637553 -0.0002640930 -0.0002643916 -0.0002646555 -0.0002648888 -0.0002650950 -0.0002652773 -0.0002654385 -0.0002655809 -0.0002657069 -0.0002658182 -0.0002659166 -0.0002660036 -0.0002660805 -0.0002661485 -0.0002662086 -0.0002662618 -0.0002663087 -0.0002663500 -0.0002663869 -0.0002664194 -0.0002664481'
  []
[]

[AuxVariables]
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [tot_force]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [tot_force]
    type = ParsedAux
    coupled_variables = 'stress_yy porepressure'
    execute_on = timestep_end
    variable = tot_force
    expression = '-stress_yy+0.6*porepressure'
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
  [grad_stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 2
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
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 0.6
    component = 2
    variable = disp_z
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

[Postprocessors]
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
  exodus = true
[]