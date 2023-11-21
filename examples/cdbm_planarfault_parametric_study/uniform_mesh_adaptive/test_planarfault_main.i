##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.1
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677
# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.1) = 0.2
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.1) * 120e6) = 185m
# Use mesh size = 25m
##########################################################

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../../../meshgenerator/cdbm/planarfault2_uniform_adaptive/planarfault2_uniform_adaptive.msh'
    []
    [./new_block_1]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y<0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block_1
        split_interface = true
        block_pairs = '0 1'
    []
    [./sidesets]
      input = split
      type = SideSetsFromNormalsGenerator
      normals = '-1 0 0
                  1 0 0
                  0 -1 0
                  0 1 0'
      new_boundary = 'left right bottom top'
    []
  []
  
  [GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y'
    
    #damping ratio
    q = 0.2
    
    #characteristic length (m)
    Dc = 0.4
  
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
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
  
    ##Compute gamma_damaged_r, xi_1
    #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
    #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
    #check struct_param.m 
  
    #coefficient of damage solid modulus
    gamma_damaged_r = 3.7150e10
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.8248
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.75
    a0 = 7.4289e9
    a1 = -2.214e10
    a2 = 2.0929e10
    a3 = -6.0672e9
  
    #diffusion coefficient #for structural stress coupling
    D = 0
    
  []
  
  [AuxVariables]
    [./resid_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_slipweakening_x]
      order = FIRST
      family = LAGRANGE
    []
    [./accel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_d]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_shear_stress]
        order = FIRST
        family = LAGRANGE
    []
    [./tangent_jump_rate]
        order = CONSTANT
        family = MONOMIAL
    []
    [./tria_area_aux]
        order = CONSTANT
        family = MONOMIAL
    []
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    #obtain parameters from MaterialRealAux, pass parameters to subApp
    # [./alpha_old]
    #   order = FIRST
    #   family = LAGRANGE
    # []
    [./B_old]
      order = FIRST
      family = LAGRANGE
    []
    [./xi_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./I2_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./lambda_old]
        order = CONSTANT
        family = MONOMIAL
    []
    [./gamma_old]
        order = CONSTANT
        family = MONOMIAL
    []
    #updated alpha, B
    [./alpha_in]
      order = CONSTANT
      family = MONOMIAL
    []
    [./B_in]
      order = CONSTANT
      family = MONOMIAL
    []
    #high-order-dummy
    [./alpha_in_dummy]
      order = FIRST
      family = MONOMIAL
    []
    [./B_in_dummy]
      order = FIRST
      family = MONOMIAL
    []
    #grad_alpha
    [./alpha_grad_x]
      order = CONSTANT
      family = MONOMIAL
    []
    [./alpha_grad_y]
      order = CONSTANT
      family = MONOMIAL
    []
    #mechanical strain rate
    [./mechanical_strain_rate]
      order = CONSTANT
      family = MONOMIAL
    []
    #track Cd
    [./track_Cd]
      order = CONSTANT
      family = MONOMIAL
    []
    #eqv plastic strain
    [./eqv_plastic_strain]
      order = CONSTANT
      family = MONOMIAL
    [] 
    #plastic work
    [./plastic_work]
      order = CONSTANT
      family = MONOMIAL
    []
  []
  
  [Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
      boundary = 'Block0_Block1'
      strain = SMALL
      generate_output='tangent_jump traction_x traction_y jump_x jump_y'
    [../]
  []
  
  
  [Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
          extra_vector_tags = 'restore_tag'
        [../]
      [../]
    [../]
  []
  
  [Problem]
    extra_tag_vectors = 'restore_tag'
  []
  
  [AuxKernels]
    [Displacment_x]
      type = ProjectionAux
      variable = disp_slipweakening_x
      v = disp_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_y]
      type = ProjectionAux
      variable = disp_slipweakening_y
      v = disp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_x]
        type = CompVarRate
        variable = vel_slipweakening_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_slipweakening_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_x]
      type = ProjectionAux
      variable = resid_slipweakening_x
      v = resid_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
      type = ProjectionAux
      variable = resid_slipweakening_y
      v = resid_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [restore_x]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_x'
      variable = 'resid_x'
      execute_on = 'TIMESTEP_END'
    []
    [restore_y]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_y'
      variable = 'resid_y'
      execute_on = 'TIMESTEP_END'
    []
    [StaticFricCoeff]
      type = FunctionAux
      variable = mu_s
      function = func_static_friction_coeff_mus
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [DynamicFricCoeff]
        type = FunctionAux
        variable = mu_d
        function = func_dynamic_friction_coeff_mud
        execute_on = 'INITIAL TIMESTEP_BEGIN'
      []
    [StrikeShearStress]
      type = FunctionAux
      variable = ini_shear_stress
      function = func_initial_strike_shear_stress
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [TJump_rate]
      type = FDCompVarRate
      variable = tangent_jump_rate
      coupled = tangent_jump
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Accel_x]
      type = FDCompVarRate2
      variable = accel_slipweakening_x
      coupled = vel_slipweakening_x
      execute_on = 'TIMESTEP_END'
    []
    [Accel_y]
      type = FDCompVarRate2
      variable = accel_slipweakening_y
      coupled = vel_slipweakening_y
      execute_on = 'TIMESTEP_END'
    []
    #obtain parameters from MaterialRealAux
    # [get_alpha_old]
    #     type = MaterialRealAux
    #     property = alpha_damagedvar
    #     variable = alpha_old
    #     execute_on = 'INITIAL TIMESTEP_BEGIN'
    # []
    # [get_B_old]
    #     type = MaterialRealAux
    #     property = B
    #     variable = B_old
    #     execute_on = 'INITIAL TIMESTEP_BEGIN'
    # []
    [get_xi_old]
        type = MaterialRealAux
        property = xi
        variable = xi_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_I2_old]
        type = MaterialRealAux
        property = I2
        variable = I2_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_mu_old]
        type = MaterialRealAux
        property = shear_modulus
        variable = mu_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_lambda_old]
        type = MaterialRealAux
        property = lambda
        variable = lambda_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [get_gamma_old]
        type = MaterialRealAux
        property = gamma_damaged
        variable = gamma_old
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #define shear strain material property (elastic) inside damage stress 
    #and compute its rate using "MaterialRateRealAux"
    [get_shear_strain_rate]
        type = MaterialRateRealAux
        property = principal_strain
        variable = mechanical_strain_rate
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #get eqv plastic strain
    [eqv_plastic_strain]
        type = MaterialRealAux
        property = eqv_plastic_strain
        variable = eqv_plastic_strain
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #get plastic work
    [plastic_work]
        type = MaterialRealAux
        property = plastic_work
        variable = plastic_work
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #fault length
    [fault_len]
        type = ConstantAux
        variable = nodal_area
        value = 25
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  []
  
  [Kernels]
    [./inertia_x]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_x
    []
    [./inertia_y]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_y
    []
    [./Reactionx]
      type = StiffPropDamping
      variable = 'disp_x'
      component = '0'
    []
    [./Reactiony]
      type = StiffPropDamping
      variable = 'disp_y'
      component = '1'
    []
  []
  
  [Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeDamageBreakageStressv2
        option = 1
        alpha_in = alpha_in_dummy
        B_in = B_in_dummy
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        output_properties = 'eps_p eps_e eps_total I1'
        outputs = exodus
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    #SlipWeakeningMultifaults ONLY supports TRIA currently!
    [./czm_mat]
        type = SlipWeakeningMultifaults
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        nodal_area = nodal_area
        mu_d = mu_d
        mu_s = mu_s
        tria_area = tria_area_aux
        boundary = 'Block0_Block1'
    [../]
    [./static_initial_stress_tensor_slipweakening]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor_slipweakening
        tensor_functions = 'func_initial_stress_xx                func_initial_strike_shear_stress      func_initial_stress_00 
                            func_initial_strike_shear_stress      func_initial_stress_yy                func_initial_stress_00
                            func_initial_stress_00                func_initial_stress_00                func_initial_stress_00'
    [../]
    [./static_initial_stress_tensor]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_initial_stress_xx             func_initial_stress_xy_const        func_initial_stress_00 
                            func_initial_stress_xy_const       func_initial_stress_yy              func_initial_stress_00
                            func_initial_stress_00             func_initial_stress_00              func_initial_stress_00'
    [../]
  []
  
  [Functions]
    #mus constant value: 0.7
    [func_static_friction_coeff_mus]
        type = ConstantFunction
        value = 0.677
    []
    #mud constant value: 0.4
    [func_dynamic_friction_coeff_mud]
        type = ConstantFunction
        value = 0.1
    []
    #Note:restrict stress variation along the fault only
    #this function is used in czm only
    [func_initial_strike_shear_stress]
      type = InitialStressXYcontmfbfs
      # type = ConstantFunction
      # value = 70e6
    []
    #this function is used in medimum
    [func_initial_stress_xy_const]
      type = ConstantFunction
      value = 70e6
    []
    [./func_initial_stress_00]
      type = ConstantFunction
      value = 0.0
    []
    [./func_initial_stress_yy]
      type = ConstantFunction
      value = -120e6
    []
    #In problems with inelasticity, the sigma11 is important
    #This is different from pure tpv205 
    [./func_initial_stress_xx]
      type = ConstantFunction
      value = -135e6
    []
  []
  
  [UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
  []
  
  [Executioner]
    type = Transient
    dt = 5e-4
    end_time = 3.0
    # num_steps = 10
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
  []
  
  #for cluster run
  [Outputs]
    exodus = true
    interval = 200
    [sample_snapshots]
      type = Exodus
      interval = 2000
    []
    [snapshots]
      type = Exodus
      interval = 2000
      overwrite = true
    []
    [checkpoints]
      type = Checkpoint
      interval = 2000
      num_files = 2
    []
  []
  
  [BCs]
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 6000
        shear_wave_speed = 3464
        boundary = right
    []
  []
  
  [MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'test_planarfault_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
  []
  
  [Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub track_Cd alpha_checked_dummy B_checked_dummy'
        variable = 'alpha_in B_in alpha_grad_x alpha_grad_y track_Cd alpha_in_dummy B_in_dummy'
        execute_on = 'TIMESTEP_BEGIN'
    []
    #we actually don't need to pass alpha and B
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate alpha_in_dummy B_in_dummy'
        variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate_sub alpha_old_dummy B_old_dummy'
        execute_on = 'TIMESTEP_BEGIN'
    []
  []