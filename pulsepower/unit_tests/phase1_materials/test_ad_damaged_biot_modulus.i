# Unit test for ElkADPorousFlowDamagedBiotModulus
# Tests that the damaged Biot modulus is computed correctly:
#   1/M = phi(d)/Kf + (alpha(d) - phi(d))/Ks
#
# Uses damaged properties from:
#   - ElkADPorousFlowDamagedBiotCoefficient: alpha(d) = 1 - g(d)*K0/Ks
#   - ElkADPorousFlowDamagedPorosity: phi(d) = phi_0 + (1-phi_0)*(1-g(d))
#
# Test parameters:
#   K0 = Ks = 50e9 Pa (so alpha = 1 - g for simplicity)
#   Kf = 2.24e9 Pa
#   phi_0 = 0.008
#
# Test cases:
#   d = 0.0: alpha = 0, phi = 0.008, 1/M = 0.008/2.24e9 + (0-0.008)/50e9 = 3.57e-12 - 1.6e-13
#   d = 0.5: alpha = 0.75, phi = 0.752, 1/M = 0.752/2.24e9 + (0.75-0.752)/50e9
#   d = 1.0: alpha ≈ 1, phi ≈ 1, 1/M = 1/2.24e9 + (1-1)/50e9 = 4.46e-10

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [biot_modulus_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [biot_coeff_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  [damage_ic]
    type = FunctionIC
    variable = d
    function = 'x'
  []
[]

[AuxKernels]
  [get_biot_modulus]
    type = ADMaterialRealAux
    variable = biot_modulus_aux
    property = PorousFlow_constant_biot_modulus_qp
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient_damaged
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_porosity]
    type = ADMaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Materials]
  # First compute damaged Biot coefficient
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = 2e-11  # 1/K0 = 1/50e9
    grain_bulk_modulus = 50e9      # Ks
    minimum_degradation = 1e-8
  []
  # Then compute damaged porosity
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
  []
  # Finally compute damaged Biot modulus using the damaged properties
  [damaged_biot_modulus]
    type = ElkADPorousFlowDamagedBiotModulus
    use_damaged_biot = true
    use_damaged_porosity = true
    fluid_bulk_modulus = 2.24e9
    grain_bulk_modulus = 50e9
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = dummy
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = dummy
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = dummy
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Postprocessors]
  [biot_modulus_left]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.1667 0.05 0'
  []
  [biot_modulus_mid]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.5 0.05 0'
  []
  [biot_modulus_right]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.8333 0.05 0'
  []
  [biot_coeff_mid]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.5 0.05 0'
  []
  [porosity_mid]
    type = PointValue
    variable = porosity_aux
    point = '0.5 0.05 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
