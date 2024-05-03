[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 400
        ny = 60
        xmin = -20000
        xmax = 20000
        ymin = -3000
        ymax = 3000
        elem_type = TRI3
    [../]
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    #add "Block0_Block1" and "Block1_Block0" interfaces
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
        add_interface_on_two_sides = true
    []
[]
  
  [GlobalParams]
  
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 3.204e10
      
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 3.204e10
  
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
    xi_1 = 0.8248
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    #coefficients
    # chi = 0.75
    a0 = 7.4289e9
    a1 = -2.214e10
    a2 = 2.0929e10
    a3 = -6.0672e9
  
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
  []
  
  #by default, subApp is using the same time step as mainApp
  [Executioner]
    type = Transient
    [TimeIntegrator]
      type = ExplicitSSPRungeKutta
      order = 3
    []
  []