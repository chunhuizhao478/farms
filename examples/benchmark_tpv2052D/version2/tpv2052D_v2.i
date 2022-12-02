# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #

# [Version 2] #
# This file serves to: #
# 1. Release the zero normal traction perturbation restriction by passing reaction forces plus/mins surfaces #
# 2. Generalize the computation of nodal triburary area along the fault surface using NodalArea function in Contact Module #

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 2000
        ny = 2000
        xmin = -50000
        xmax = 50000
        ymin = -50000
        ymax = 50000
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y>0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
    mu_d = 0.525
[]

[AuxVariables]
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_x]
    []
    [./vel_y]
    []
    [./accel_y]
    []
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_primary_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_primary_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_secondary_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_secondary_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y normal_traction tangent_traction normal_jump tangent_jump'
        save_in_master = 'resid_primary_x resid_primary_y'
        save_in_slave  = 'resid_secondary_x  resid_secondary_y' 
    [../]
[]


[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            add_variables = true
            planar_formulation = PLANE_STRAIN
            generate_output = 'stress_xx stress_yy stress_xy'
        [../]
        [../]
    [../]
[]

[AuxKernels]
    [velocity_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
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
        type = SlipWeakeningFriction2dv2
        boundary = 'Block0_Block1'
        react_x = resid_primary_x
        react_y = resid_primary_y
        react_neighbor_x = resid_secondary_x
        react_neighbor_y = resid_secondary_y
        nodal_area = nodal_area
    [../]
[]

[UserObjects]
    [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = Block0_Block1
      execute_on = 'initial TIMESTEP_BEGIN'
    [../]
[]

[Executioner]
    type = Transient
    dt = 0.005
    end_time = 2
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 10
[]