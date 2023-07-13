#This is the main file for solving for rate-and-state friction

#The following main quantities are defined as material property and are declared in "CZMComputeLocalTractionBaseRSF2D":

#CZMComputeLocalTractionBaseRSF2D <- CZMComputeLocalTractionTotalBaseRSF2D <- RateStateFrictionLaw2Dv7

#traction_strike (traction along strike direction)
#sliprate_strike (sliprate along strike direction)
#slip_strike     (slip along strike direction)
#statevar        (statevar)

#The information stores in the blocks attached to the primary surface

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 300
        ny = 300
        xmin = -15000
        xmax = 15000
        ymin = -15000
        ymax = 15000
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
    ##primary variable
    displacements = 'disp_x disp_y'
    
    ##damping ratio 
    q = 0.1
    
    ##element length (m)
    len = 100
    
    ##rate-and-state coefficients
    f_o = 0.6
    rsf_a = 0.008
    rsf_b = 0.012
    rsf_L = 0.02
    delta_o = 1e-6

    ##initial normal traction (Pa)
    Tn_o = 120e6

    ##initial shear traction (Pa)
    Ts_o = 75e6

    ##initial sliprate (m/s)
    Vini = 1e-12

    ##initial state variable
    statevarini = 1.606238999213454e9
[]

[Variables]
    [./disp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    ##restoration force (pass into material kernel)
    #^the vector takes sub app evaluation values
    [./resid_rsf_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_rsf_y]
        order = FIRST
        family = LAGRANGE
    [../]
    
    ##velocity vector (output)
    #^use 1st-order finite difference to approx the derivative ("CoupledDot() is working only with nonlinear variable")
    #^^This is time derivative of total displacement
    [./vel_rsf_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_rsf_y]
        order = FIRST
        family = LAGRANGE
    []

    ##initial shear stress
    #^take aux variable from "FunctionAux"
    [./ini_shear_stress_perturb]
        order = FIRST
        family = LAGRANGE
    []

    #tangent jump
    #^take aux variable from CZM outputs
    [./tangent_jump]
        order = CONSTANT
        family = MONOMIAL
    []

    #tangent jump rate
    #^use 1st-order finite difference to approx the derivative
    [./tangent_jump_rate]
        order = CONSTANT
        family = MONOMIAL
    []

    #slip rate
    #^take aux variable from sub App for output/visualization purpose
    [./sliprate_strike]
        order = CONSTANT
        family = MONOMIAL
    []

    #background displacement
    #^constant velocity along x direction
    [./vel_const_x]
        order = FIRST
        family = LAGRANGE
    []

    #^constant velocity along y direction (zero, not used)
    [./vel_const_y]
        order = FIRST
        family = LAGRANGE
    []

    ##total displacement
    #^displacement vector takes disp + background disp
    [./disp_rsf_x]
        order = FIRST
        family = LAGRANGE
    [] 
    [./disp_rsf_y]
        order = FIRST
        family = LAGRANGE
    [] 

[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output = 'tangent_jump'
    [../]
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          planar_formulation = PLANE_STRAIN
          displacements = 'disp_x disp_y'
        [../]
      [../]
    [../]
[]

[AuxKernels]
    ##Compute total displacement (disp_rsf_x = disp_x + vel_const * _t) 
    #^disp_rsf needs to be updated at each time step end
    [Disp_Total_x]
        type = CompDispRateState
        variable = disp_rsf_x
        current_disp = disp_x
        vel_const = vel_const_x
        execute_on = 'TIMESTEP_END'
    []
    [Disp_Total_y]
        type = CompDispRateState
        variable = disp_rsf_y
        current_disp = disp_y
        vel_const = vel_const_y
        execute_on = 'TIMESTEP_END'
    []

    ##Compute velocity (output) 
    #^Compute the time derivative of total displacement ???
    [Vel_x]
        type = FDCompVarRatev2
        variable = vel_rsf_x
        coupled = disp_rsf_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = FDCompVarRatev2
        variable = vel_rsf_y
        coupled = disp_rsf_y
        execute_on = 'TIMESTEP_END'
    []

    ##Add background tangent jump
    [TJump_addbd]
        type = CompTJumpRateState
        variable = tangent_jump
        coupled = tangent_jump
        sliprate_bd = 1e-12
        execute_on = 'TIMESTEP_END'
    []

    ##Compute tangent jump rate
    [TJump_rate]
        type = FDCompVarRatev2
        variable = tangent_jump_rate
        coupled = tangent_jump
        execute_on = 'TIMESTEP_END'
    []
    
    ##initial shear traction
    [StrikeShearStress]
        type = FunctionAux
        variable = ini_shear_stress_perturb
        function = func_initial_strike_shear_stress
        execute_on = 'TIMESTEP_BEGIN'
    []

    ##output
    [outputsliprate]
        type = MaterialRealAux
        property = sliprate_strike
        variable = sliprate_strike
        boundary = 'Block0_Block1'
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
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    [./czm_mat]
        type = RateStateFrictionLaw2Dv7
        reaction_rsf_x = resid_rsf_x
        reaction_rsf_y = resid_rsf_y
        Ts_perturb = ini_shear_stress_perturb
        boundary = 'Block0_Block1'
    [../]
[]

#^For testing, we close the stress perturbation
[Functions]
    [./func_initial_strike_shear_stress]
        # type = InitialStrikeShearStressPerturbRSF2D
        type = ConstantFunction
        value = 0
    []
[]

[Executioner]
    type = Transient
    dt = 0.0025
    end_time = 3.0
    num_steps = 5
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]

[MultiApps]
    #^This subApp here is only used for residual evaluation
    [./sub_app0]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'ratestate2d_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    [pull_resid0]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app0
        source_variable = 'resid_sub_x resid_sub_y'
        variable = 'resid_rsf_x resid_rsf_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp0]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app0
        source_variable = 'disp_rsf_x disp_rsf_y'
        variable = 'disp_sub_x disp_sub_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

# [Debug]
#     show_execution_order = ALWAYS
# []

[ICs]
    [./initial_vel_total_x_primary]
        type = ConstantIC
        variable = vel_const_x
        value = 5e-13
        block = 0
    []
    [./initial_vel_total_x_secondary]
        type = ConstantIC
        variable = vel_const_x
        value = -5e-13
        block = 1
    []
[]