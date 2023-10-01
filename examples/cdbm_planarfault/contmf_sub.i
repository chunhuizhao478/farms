#sample test geometry
[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 800
    ny = 200
    xmin = -20000
    xmax = 20000
    ymin = -5000
    ymax = 5000
    elem_type = TRI3
  []
  [./new_block]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'y<0'
    block_id = 1
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
    split_interface = true
  []
[]
  
[GlobalParams]

    ##----continuum damage breakage model----##

    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8

    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.9

    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.5

    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.5

    #<coefficient gives positive damage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_d = 40e5

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_B = 40e7

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    C_BH = 0

    #<coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4
    C_1 = 0

    #<coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3

    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.8248

    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m

    #Note: "computeAlphaCr" needs to change every time the related parameters changed

    #coefficients
    a0 = 4.9526e9
    a1 = -1.8888e10
    a2 = 2.3960e10
    a3 = -1.0112e10

    #diffusion coefficient #self-defined value
    D = 100
    
  []

  [Variables]
    [alpha_sub]
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
    #
    [B_sub]
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
        type = DamageVarForcingFunc
        alpha_old = alpha_old
        B_old = B_old
        xi_old = xi_old
        I2_old = I2_old
        variable = alpha_sub
    []
  []

  [AuxKernels]
    [compute_B]
      type = BreakageVarUpdate
      variable = B_sub
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      execute_on = 'TIMESTEP_BEGIN'
    []
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
  []

  #by default, subApp is using the same time step as mainApp
  [Executioner]
    type = Transient
    [TimeIntegrator]
      type = ExplicitSSPRungeKutta
      order = 3
    []
  []

  # [Outputs]
  #   exodus = true
  #   interval = 200
  # []