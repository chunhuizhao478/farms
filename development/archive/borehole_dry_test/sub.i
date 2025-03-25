[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  './meshfile/borehole_wofaults.msh'
  []
[]
  
[GlobalParams]
  
  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  # lambda_o = 14e9
    
  #initial shear modulus value (second lame constant) [Pa]
  # shear_modulus_o = 14e9

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.929

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.929

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
  Cd_constant = 1e4

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
  xi_1 = 0.7631

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  # #coefficients
  # chi = 0.75
  a0 = 2.788e9
  a1 = -8.6937e9
  a2 = 8.4215e9
  a3 = -2.3806e9

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
    #grad_alpha
    [alpha_grad_x_sub]
      order = CONSTANT
      family = MONOMIAL
    []
    [alpha_grad_y_sub]
      order = CONSTANT
      family = MONOMIAL
    []
[]
  
[Kernels]
    [./timederivative_alpha]
        type = ADTimeDerivative
        variable = alpha_sub
    []
    [./alpha_forcing_func]
        type = ADDamageVarForcingFunc
        alpha_old = alpha_old
        B_old = B_old
        xi_old = xi_old
        I2_old = I2_old
        variable = alpha_sub
    []
    [./timederivative_B]
      type = ADTimeDerivative
      variable = B_sub
    []
    [./B_forcing_func]
        type = ADBreakageVarForcingFunc
        variable = B_sub
        alpha_old = alpha_old
        B_old = B_old
        xi_old = xi_old
        I2_old = I2_old
        mu_old = mu_old
        gamma_old = gamma_old
        lambda_old = lambda_old
    []
[]

[AuxKernels]
    [check_alpha]
        type = CheckAlpha
        coupled = alpha_sub
        variable = alpha_checked
        execute_on = 'TIMESTEP_END'
    []
    [check_B]
        type = CheckB
        coupled = B_sub
        variable = B_checked
        execute_on = 'TIMESTEP_END'
    []
[]
  
#by default, subApp is using the same time step as mainApp
[Executioner]
    type = Transient
    [TimeIntegrator]
      type = ImplicitEuler
    []
[]