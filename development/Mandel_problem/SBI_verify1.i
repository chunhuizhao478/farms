

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 500
  ny = 250
  xmin = -50
  xmax = 50
  ymin = 0
  ymax = 25
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]


[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

# ========================================================================
# BOUNDARY CONDITIONS - Apply these functions to your boundaries
# ========================================================================
[BCs]
  # Apply horizontal displacement on bottom boundary (z=0)
  [ux_bottom]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'  # change to match your mesh boundary name
    function = ux_bc
  []

  # Apply vertical displacement on bottom boundary (z=0)
  [uz_bottom]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = uz_bc
  []

[]



[Functions]
  # ========================================================================
  # HORIZONTAL DISPLACEMENT ux(x,t) = U0 * phi(x) * g(t)
  # ux(x,t) = U0 * exp(-((x-x0)^2)/Lx^2) * (t/t0)*exp(-t/t0)
  # ========================================================================
  [ux_bc]
    type = ParsedFunction
    expression = 'if(t <= 0, 0, U0 * exp(-((x - x0)^2) / (Lx^2)) * (t/t0) * exp(-t/t0))'
    symbol_names = 'U0   x0   Lx   t0'
    symbol_values = '5  0.0  2.0  1.0'
    # U0: displacement amplitude (5*U0 from MATLAB)
    # x0: center of Gaussian patch
    # Lx: half-width of localization
    # t0: time scale of pulse (set to ~t_end/5)
  []

  # ========================================================================
  # VERTICAL DISPLACEMENT uz(x,t) = beta*U0 * phi(x) * g(t)
  # uz(x,t) = beta*U0 * exp(-((x-x0)^2)/Lx^2) * (t/t0)*exp(-t/t0)
  # ========================================================================
  [uz_bc]
    type = ParsedFunction
    expression = 'if(t <= 0, 0, U0 * exp(-((x - x0)^2) / (Lx^2)) * (t/t0) * exp(-t/t0))'
    symbol_names = 'U0   x0   Lx   t0'
    symbol_values = '1  0.0  2.0  1.0'
    # U0: displacement amplitude (beta*U0 = 0.5*1.0 from MATLAB)
    # x0: center of Gaussian patch
    # Lx: half-width of localization
    # t0: time scale of pulse
  []
[]



[AuxVariables]
  [stress_xx]
    order = FIRST
    family = MONOMIAL
  []
  [stress_yy]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 1
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
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
[]

[Materials]
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '0.428571 1'
    # bulk modulus is lambda + 2*mu/3 = 0.5 + 2*0.75/3 = 1
    fill_method = symmetric_isotropic
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
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
  end_time = 5
  dt = 0.025
[]

[Outputs]
  exodus = true
[]