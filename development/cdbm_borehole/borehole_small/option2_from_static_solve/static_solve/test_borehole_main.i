##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4 [Check!]
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677

# under the stress state, 
# -135MPa 70MPa; 70MPa -120e6
# the principal stress is 198e6 

# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.5
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.4) * 198e6) = 233.6m

# Use mesh size = 25m, resolved by more than 7 elements [Check!]

# Use 1% overstress

# Use dt = 2e-4
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../../../meshgenerator/cdbm/borehole_wofaults_small/borehole_wofaults.msh'
  []
  [./subdomain_id] 
    input = msh
    type = ElementSubdomainIDGenerator 
    element_ids = '
    21 265 499 955 995 1227 1543 2001 2397 2413 2747 2827 3071 3855 4137 4305 4315 4497 4521 4575 4693 4787 5059 5633 5867 6125 6219 6403 6487 6831 6929 7319 7365 7719 7833 7849 8029 8043 8071 8299 8305 8595 8596 9345 9375 9423 9603 9667 9939 9943 10231 10535 10715 10729 10905 11069 11221 11289 11401 11413 11523 11619 11853 12893 12967 13047 13063 13113 13199 13327 13353 13541 13643 13711 13721 13737 13769 13839 13939 13955 20 264 498 954 994 1226 1542 2000 2396 2412 2746 2826 3070 3854 4136 4304 4314 4496 4520 4574 4692 4786 5058 5632 5866 6124 6218 6402 6486 6830 6928 7318 7364 7718 7832 7848 8028 8042 8070 8298 8304 8594 8597 9344 9374 9422 9602 9666 9938 9942 10230 10534 10714 10728 10904 11068 11220 11288 11400 11412 11522 11618 11852 12892 12966 13046 13062 13112 13198 13326 13352 13540 13642 13710 13720 13736 13768 13838 13938 13954 23 267 501 953 993 1225 1541 2003 2395 2411 2745 2825 3073 3853 4135 4307 4317 4495 4519 4577 4691 4785 5057 5635 5869 6127 6221 6401 6489 6829 6927 7321 7363 7717 7831 7851 8027 8041 8073 8297 8303 8593 8599 9343 9377 9425 9601 9665 9937 9945 10229 10533 10713 10731 10903 11071 11223 11287 11403 11415 11525 11617 11855 12891 12969 13049 13065 13115 13201 13325 13355 13539 13641 13713 13719 13739 13767 13837 13941 13953 22 266 500 952 992 1224 1540 2002 2394 2410 2744 2824 3072 3852 4134 4306 4316 4494 4518 4576 4690 4784 5056 5634 5868 6126 6220 6400 6488 6828 6926 7320 7362 7716 7830 7850 8026 8040 8072 8296 8302 8592 8598 9342 9376 9424 9600 9664 9936 9944 10228 10532 10712 10730 10902 11070 11222 11286 11402 11414 11524 11616 11854 12890 12968 13048 13064 13114 13200 13324 13354 13538 13640 13712 13718 13738 13766 13836 13940 13952 812 813 1414 1415 1676 1677 1684 1685 1720 1721 1838 1839 2072 2073 2086 2087 2196 2197 2240 2241 2594 2595 2774 2775 3724 3725 3818 3819 3892 3893 4002 4003 4412 4413 4616 4617 4840 4841 4888 4889 4930 4931 5088 5089 5148 5149 5328 5329 5636 5637 5670 5671 5772 5773 5984 5985 6228 6229 6358 6359 6392 6393 6560 6561 6590 6591 6738 6739 6770 6771 7098 7099 7130 7131 7434 7435 7548 7549 8010 8011 8022 8023 8162 8163 8406 8407 8412 8413 8538 8539 8576 8577 8938 8939 9036 9037 9354 9355 9468 9469 9694 9695 10090 10091 10686 10687 10734 10735 11074 11075 11580 11581 11606 11607 12616 12617 13120 13121 13280 13281 13302 13303 13422 13423 13486 13487 13840 13841 14124 14125 14258 14259 14276 14277 14370 14371 14388 14389 15070 15071 
    '
    subdomain_ids = '
    1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 
    '
  []
  [./sidesets_leftbc]
    input = subdomain_id
    type = SideSetsAroundSubdomainGenerator
    block = 1000
    new_boundary = 'left'
    external_only = true
  []
  [./sidesets_rightbc]
    input = sidesets_leftbc
    type = SideSetsAroundSubdomainGenerator
    block = 1001
    new_boundary = 'right'
    external_only = true
  []
  [./sidesets_topbc]
    input = sidesets_rightbc
    type = SideSetsAroundSubdomainGenerator
    block = 1002
    new_boundary = 'top'
    external_only = true
  []
  [./sidesets_bottombc]
    input = sidesets_topbc
    type = SideSetsAroundSubdomainGenerator
    block = 1003
    new_boundary = 'bottom'
    external_only = true
  []
  [./sidesets_boreholebc]
    input = sidesets_bottombc
    type = SideSetsAroundSubdomainGenerator
    block = 1004
    new_boundary = 'borehole'
    external_only = true
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0  -1000  0'
    new_boundary = corner_ptr1
    input = sidesets_boreholebc
  []
  [./extranodeset2]
      type = ExtraNodesetGenerator
      coord = '1000  0  0'
      new_boundary = corner_ptr2
      input = extranodeset1
  []
  [./extranodeset3]
      type = ExtraNodesetGenerator
      coord = '0   1000  0'
      new_boundary = corner_ptr3
      input = extranodeset2
  []
  [./extranodeset4]
      type = ExtraNodesetGenerator
      coord = '-1000  0  0'
      new_boundary = corner_ptr4
      input = extranodeset3
  []
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  # #damping ratio
  # q = 0.2
  
  # #characteristic length (m)
  # # Dc = 0.4

  # ##----continuum damage breakage model----##
  # #initial lambda value (first lame constant) [Pa]
  # lambda_o = 3.204e10
    
  # #initial shear modulus value (second lame constant) [Pa]
  # shear_modulus_o = 3.204e10

  # #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  # xi_0 = -1.0

  # #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  # xi_d = -1.1

  # #<strain invariants ratio: maximum allowable value>: set boundary
  # #Xu_etal_P15-2D
  # #may need a bit space, use 1.5 as boundary
  # xi_max = 1.8

  # #<strain invariants ratio: minimum allowable value>: set boundary
  # #Xu_etal_P15-2D
  # xi_min = -1.8

  # #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  # C_g = 1e-10

  # #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  # m1 = 10

  # #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  # m2 = 1

  # ##Compute gamma_damaged_r, xi_1
  # #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
  # #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
  # #check struct_param.m 

  # #coefficient of damage solid modulus
  # gamma_damaged_r = 3.49061e10

  # #critical point of three phases (strain invariants ratio vs damage)
  # xi_1 = 0.73205

  # ##Compute parameters in granular states
  # #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  # #check struct_param.m

  # #--------------------------------------------------------------------------------#
  # #Note: "computeAlphaCr" needs to change every time the related parameters changed
  # #--------------------------------------------------------------------------------#

  # # #coefficients
  # # chi = 0.75
  # a0 = 3.92486e9
  # a1 = -1.51257e10
  # a2 = 1.93525e10
  # a3 = -8.21570e9

  # #diffusion coefficient #for structural stress coupling
  # D = 0

  # #pressure analytical solution
  # #Reference: Injection-induced seismicity: Poroelastic and earthquake nucleation effects (P. Segall1 and S. Lu2)
  # # effec_sts_coeff = 0.31 #
  # # flux_q = 5e2 #kg/s
  # # density_rho_0 = 1e3 #kg/m^3
  # # permeability_k = 3e-12 #m^2
  # # viscosity_eta = 0.4e-3 #Pa s
  # biotcoeff_alpha = 0.31 #-
  # # undrained_nu_u = 0.3  #-
  # # shear_modulus_mu = 32.04e9 #Pa
  # # drained_nu = 0.25 #-
  
[]

[AuxVariables]
  [./stress_xx_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./stress_yy_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./strain_yy_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xi_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
#   [./resid_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./resid_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   [./resid_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   [./disp_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./disp_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./vel_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./vel_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./mu_s]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./mu_d]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./ini_shear_stress]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./tangent_jump_rate]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./tria_area_aux]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./nodal_area]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   #obtain parameters from MaterialRealAux, pass parameters to subApp
#   # [./alpha_old]
#   #   order = FIRST
#   #   family = LAGRANGE
#   # []
#   [./B_old]
#     order = FIRST
#     family = LAGRANGE
#   []
#   [./xi_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./I2_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./mu_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./lambda_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./gamma_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   #updated alpha, B
#   [./alpha_in]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./B_in]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #high-order-dummy
#   [./alpha_in_dummy]
#     order = FIRST
#     family = MONOMIAL
#   []
#   [./B_in_dummy]
#     order = FIRST
#     family = MONOMIAL
#   []
#   #grad_alpha
#   [./alpha_grad_x]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./alpha_grad_y]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #mechanical strain rate
#   [./mechanical_strain_rate]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #track Cd
#   [./track_Cd]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #pressure
#   [pressure]
#     order = FIRST
#     family = LAGRANGE
#   []
# []

# [Modules/TensorMechanics/CohesiveZoneMaster]
#   [./czm_ik]
#     boundary = 'czm'
#     strain = SMALL
#     generate_output='tangent_jump normal_jump'
#   [../]
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'strain_xx strain_xy strain_yy'
        # extra_vector_tags = 'restore_tag'
      [../]
    [../]
  [../]
[]

# [Problem]
#   extra_tag_vectors = 'restore_tag'
# []

[AuxKernels]
  [get_stress_xx]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 0
    variable = stress_xx_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_xy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 1
    variable = stress_xy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_yy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 1
    variable = stress_yy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  #
  [get_strain_xx]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 0
    j = 0
    variable = strain_xx_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_xy]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 0
    j = 1
    variable = strain_xy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_strain_yy]
    type = MaterialRankTwoTensorAux
    property = mechanical_strain
    i = 1
    j = 1
    variable = strain_yy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_xi]
    type = CompXi
    strainxx = strain_xx_saved
    strainxy = strain_xy_saved
    strainyy = strain_yy_saved
    variable = xi_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
#   [Residual_x]
#     type = ProjectionAux
#     variable = resid_slipweakening_x
#     v = resid_x
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Residual_y]
#     type = ProjectionAux
#     variable = resid_slipweakening_y
#     v = resid_y
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [restore_x]
#     type = TagVectorAux
#     vector_tag = 'restore_tag'
#     v = 'disp_x'
#     variable = 'resid_x'
#     execute_on = 'TIMESTEP_END'
#   []
#   [restore_y]
#     type = TagVectorAux
#     vector_tag = 'restore_tag'
#     v = 'disp_y'
#     variable = 'resid_y'
#     execute_on = 'TIMESTEP_END'
#   []
#   [StaticFricCoeff]
#     type = FunctionAux
#     variable = mu_s
#     function = func_static_friction_coeff_mus
#     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   [DynamicFricCoeff]
#       type = FunctionAux
#       variable = mu_d
#       function = func_dynamic_friction_coeff_mud
#       execute_on = 'INITIAL TIMESTEP_BEGIN'
#     []
#   [StrikeShearStress]
#     type = FunctionAux
#     variable = ini_shear_stress
#     function = func_initial_strike_shear_stress
#     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   # [TJump_rate]
#   #   type = FDCompVarRate
#   #   variable = tangent_jump_rate
#   #   coupled = tangent_jump
#   #   execute_on = 'TIMESTEP_BEGIN'
#   # []
#   #obtain parameters from MaterialRealAux
#   # [get_alpha_old]
#   #     type = MaterialRealAux
#   #     property = alpha_damagedvar
#   #     variable = alpha_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_B_old]
#   #     type = MaterialRealAux
#   #     property = B
#   #     variable = B_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_xi_old]
#   #     type = MaterialRealAux
#   #     property = xi
#   #     variable = xi_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_I2_old]
#   #     type = MaterialRealAux
#   #     property = I2
#   #     variable = I2_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_mu_old]
#   #     type = MaterialRealAux
#   #     property = shear_modulus
#   #     variable = mu_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_lambda_old]
#   #     type = MaterialRealAux
#   #     property = lambda
#   #     variable = lambda_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_gamma_old]
#   #     type = MaterialRealAux
#   #     property = gamma_damaged
#   #     variable = gamma_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # #define shear strain material property (elastic) inside damage stress 
#   # #and compute its rate using "MaterialRateRealAux"
#   # [get_shear_strain_rate]
#   #     type = MaterialRateRealAux
#   #     property = principal_strain
#   #     variable = mechanical_strain_rate
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   #fault length
#   [fault_len]
#       type = ConstantAux
#       variable = nodal_area
#       value = 25
#       execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   #initial is zero, DO NOT add initial on execute_on
#   [./aux_func_pressure]
#     type = ConstantAux
#     variable = pressure
#     value = 0
#     execute_on = 'TIMESTEP_BEGIN'
#   [../]
[]

# [Kernels]
#   [./inertia_x]
#     type = InertialForce
#     use_displaced_mesh = false
#     variable = disp_x
#   []
#   [./inertia_y]
#     type = InertialForce
#     use_displaced_mesh = false
#     variable = disp_y
#   []
#   [./Reactionx]
#     type = StiffPropDamping
#     variable = 'disp_x'
#     component = '0'
#   []
#   [./Reactiony]
#     type = StiffPropDamping
#     variable = 'disp_y'
#     component = '1'
#   []
# []

[Materials]
  # #damage breakage model
  # [stress_medium]
  #   type = ComputeDamageBreakageStressv3pressurev2
  #   alpha_in = alpha_in_dummy
  #   B_in = B_in_dummy
  #   alpha_grad_x = alpha_grad_x
  #   alpha_grad_y = alpha_grad_y
  #   pressure = pressure
  #   output_properties = 'eps_p eps_e eps_total I1'
  #   outputs = exodus
  # []
  [stress_medium]
      type = ComputeLinearElasticStress
  []
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
  []
  # [density]
  #     type = GenericConstantMaterial
  #     prop_names = density
  #     prop_values = 2670
  # []
  #SlipWeakeningMultifaults ONLY supports TRIA currently!
  # [./czm_mat]
  #     type = SlipWeakeningMultifaultsPressure
  #     disp_slipweakening_x     = disp_slipweakening_x
  #     disp_slipweakening_y     = disp_slipweakening_y
  #     reaction_slipweakening_x = resid_slipweakening_x
  #     reaction_slipweakening_y = resid_slipweakening_y
  #     nodal_area = nodal_area
  #     mu_d = mu_d
  #     mu_s = mu_s
  #     tria_area = tria_area_aux
  #     pressure = pressure
  #     boundary = 'czm'
  # [../]
  # [./static_initial_stress_tensor_slipweakening]
  #     type = GenericFunctionRankTwoTensor
  #     tensor_name = static_initial_stress_tensor_slipweakening
  #     tensor_functions = 'func_initial_stress_xx         func_initial_strike_shear_stress           func_initial_stress_00 
  #     func_initial_strike_shear_stress         func_initial_stress_yy           func_initial_stress_00
  #                         func_initial_stress_00         func_initial_stress_00           func_initial_stress_00'
  # [../]
  # [./static_initial_stress_tensor]
  #     type = GenericFunctionRankTwoTensor
  #     tensor_name = static_initial_stress_tensor
  #     tensor_functions = 'func_initial_stress_xx             func_initial_stress_xy_const        func_initial_stress_00 
  #                         func_initial_stress_xy_const       func_initial_stress_yy              func_initial_stress_00
  #                         func_initial_stress_00             func_initial_stress_00              func_initial_stress_00'
  # [../]
[]

[Functions]
  # #mus constant value: 0.7
  # [func_static_friction_coeff_mus]
  #     type = ConstantFunction
  #     value = 0.677
  # []
  # #mud constant value: 0.4
  # [func_dynamic_friction_coeff_mud]
  #     type = ConstantFunction
  #     value = 0.4
  # []
  # #Note:restrict stress variation along the fault only
  # #this function is used in czm only
  # [func_initial_strike_shear_stress]
  #   # type = InitialStressXYnetwork_old
  #   type = ConstantFunction
  #   value = 0.0
  # []
  # #this function is used in medimum
  # [func_initial_stress_xy_const]
  #   type = ConstantFunction
  #   value = 0.0
  # []
  # [./func_initial_stress_00]
  #   type = ConstantFunction
  #   value = 0.0
  # []
  [./func_initial_stress_neg]
    type = ConstantFunction
    value = -100e6
  []
  #In problems with inelasticity, the sigma11 is important
  #This is different from pure tpv205 
  [./func_initial_stress_pos]
    type = ConstantFunction
    value = 100e6
  []
  # #pressure
  # [func_pressure_x]
  #   type = InitialStressXYPressureBorehole_x_fast
  # []
  # [func_pressure_y]
  #   type = InitialStressXYPressureBorehole_y_fast
  # []
[]

# [UserObjects]
#   [recompute_residual_tag]
#       type = ResidualEvaluationUserObject
#       vector_tag = 'restore_tag'
#       force_preaux = true
#       execute_on = 'TIMESTEP_END'
#   []
# []

[Preconditioning]
  [smp]
    full = true
    type = SMP
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       NONZERO               1e-15'
  solve_type = NEWTON
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 1
  # [sample_snapshots]
  #   type = Exodus
  #   interval = 2000
  # []
  # [snapshots]
  #   type = Exodus
  #   interval = 2000
  #   overwrite = true
  # []
  # [checkpoints]
  #   type = Checkpoint
  #   interval = 2000
  #   num_files = 2
  # []
[]

[BCs]
  [Pressure_top]
    type = Pressure
    boundary = top
    variable = disp_y
    factor = 100e6
  []
  [Pressure_bottom]
    type = Pressure
    boundary = bottom
    variable = disp_y
    factor = 100e6
  []
  [Pressure_left]
    type = Pressure
    boundary = left
    variable = disp_x
    factor = 100e6
  []
  [Pressure_right]
    type = Pressure
    boundary = right
    variable = disp_x
    factor = 100e6
  []
  [./fix_cptr1_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr1
    value = 0
  []
  # [./fix_cptr1_y]
  #     type = DirichletBC
  #     variable = disp_y
  #     boundary = corner_ptr1
  #     value = 0
  # []
  # [./fix_cptr2_x]
  #     type = DirichletBC
  #     variable = disp_x
  #     boundary = corner_ptr2
  #     value = 0
  # []
  [./fix_cptr2_y]
      type = DirichletBC
      variable = disp_y
      boundary = corner_ptr2
      value = 0
  []
  [./fix_cptr3_x]
      type = DirichletBC
      variable = disp_x
      boundary = corner_ptr3
      value = 0
  []
  # [./fix_cptr3_y]
  #     type = DirichletBC
  #     variable = disp_y
  #     boundary = corner_ptr3
  #     value = 0
  # []
  # [./fix_cptr4_x]
  #     type = DirichletBC
  #     variable = disp_x
  #     boundary = corner_ptr4
  #     value = 0
  # []
  [./fix_cptr4_y]
      type = DirichletBC
      variable = disp_y
      boundary = corner_ptr4
      value = 0
  []
  # [fix_top_x]
  #   type = DirichletBC
  #   boundary = top
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_left_y]
  #   type = DirichletBC
  #   boundary = left
  #   variable = disp_y
  #   value = 0
  # []
  # [fix_bottom_x]
  #   type = DirichletBC
  #   boundary = bottom
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_bottom_y]
  #   type = DirichletBC
  #   boundary = bottom
  #   variable = disp_y
  #   value = 0
  # []
  # [fix_right_x]
  #   type = DirichletBC
  #   boundary = right
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_right_y]
  #   type = DirichletBC
  #   boundary = right
  #   variable = disp_y
  #   value = 0
  # []
[]

# [MultiApps]
#   [./sub_app]
#       type = TransientMultiApp
#       positions = '0 0 0'
#       input_files = 'test_borehole_sub.i'
#       execute_on = 'TIMESTEP_BEGIN'
#   [../]
# []

# [Transfers]
#   [pull_resid]
#       type = MultiAppCopyTransfer
#       from_multi_app = sub_app
#       source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub track_Cd alpha_checked_dummy B_checked_dummy'
#       variable = 'alpha_in B_in alpha_grad_x alpha_grad_y track_Cd alpha_in_dummy B_in_dummy'
#       execute_on = 'TIMESTEP_BEGIN'
#   []
#   #we actually don't need to pass alpha and B
#   [push_disp]
#       type = MultiAppCopyTransfer
#       to_multi_app = sub_app
#       source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate alpha_in_dummy B_in_dummy'
#       variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate_sub alpha_old_dummy B_old_dummy'
#       execute_on = 'TIMESTEP_BEGIN'
#   []
# []