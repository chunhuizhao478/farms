##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.1
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677
# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.514
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.1) * 120e6) = 185m
# Use mesh size = 25m

# Diffusion Length Scale D = 5e5
# sqrt(5e5*185/3464) = 163. using 6~7, 25m mesh to resolve it

# Check CFL condition 0.1 * 25 / 6000 ~ 0.0001s
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
  
    ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 3.204e10
    
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 3.204e10

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.6

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.6

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #<coefficient gives positive damage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #under slow strain rate < low strain rate threshold
  # C_d_min = 10

  #if option 2, use Cd_constant
  Cd_constant = 1e5

  #power-law correction
  #index
  m = 0.9

  #low strain rate threshold
  # mechanical_strain_rate_threshold = -1e-4

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = 10 

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  CBCBH_multiplier = 10

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = 300

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = 0.05

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = 0.03 #1e-3

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.936229

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  # #coefficients
  # chi = 0.75
  a0 = 6.37464e+09
  a1 = -2.63965e+10
  a2 = 3.45711e+10
  a3 = -1.45789e+10

  #diffusion coefficient #for structural stress coupling
  D = 0
  
[]

[Variables]
  [alpha_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxVariables]
  [alpha_old]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_old]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #
  [xi_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [I2_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [mu_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [lambda_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [gamma_old]
      order = CONSTANT
      family = MONOMIAL
  []
  #checked
  [alpha_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  #high-order-dummy
  [alpha_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #grad_alpha
  [alpha_grad_x_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [alpha_grad_y_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #mechanical strain rate sub
  [mechanical_strain_rate_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #Track Cd
  [track_Cd]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [./timederivative_alpha]
      type = TimeDerivative
      variable = alpha_sub
  []
  [./alpha_forcing_func]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub
      healing = true
  []
  [./timederivative_B]
    type = TimeDerivative
    variable = B_sub
  []
  [./B_forcing_func]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
  #-high-order-dummy-#
  [./timederivative_alpha_dummy]
    type = TimeDerivative
    variable = alpha_sub_dummy
  []
  [./alpha_forcing_func_dummy]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub_dummy
      healing = true
  []
  [./timederivative_B_dummy]
    type = TimeDerivative
    variable = B_sub_dummy
  []
  [./B_forcing_func_dummy]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub_dummy
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
[]

[AuxKernels]
  [check_alpha]
    type = CheckAlphaB
    coupled = alpha_sub
    variable = alpha_checked
    execute_on = 'TIMESTEP_END'
  []
  [check_B]
    type = CheckAlphaB
    coupled = B_sub
    variable = B_checked
    execute_on = 'TIMESTEP_END'
  []
  #-high-order-dummy-#
  [check_alpha_dummy]
    type = CheckAlphaB
    coupled = alpha_sub_dummy
    variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  [check_B_dummy]
    type = CheckAlphaB
    coupled = B_sub_dummy
    variable = B_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  #
  [alpha_grad_x_sub]
    type = VariableGradientComponent
    variable = alpha_grad_x_sub
    component = x
    gradient_variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
 []
 [alpha_grad_y_sub]
   type = VariableGradientComponent
   variable = alpha_grad_y_sub
   component = y
   gradient_variable = alpha_checked_dummy
   execute_on = 'TIMESTEP_END'
  []
[]

#by default, subApp is using the same time step as mainApp
[Executioner]
  type = Transient
  [TimeIntegrator]
    type = ExplicitSSPRungeKutta
    order = 3
  []
[]