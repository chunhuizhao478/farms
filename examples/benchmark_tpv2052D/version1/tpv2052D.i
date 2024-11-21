# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 100
        ny = 100
        xmin = -50000
        xmax = 50000
        ymin = -50000
        ymax = 50000
        elem_type = QUAD9
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

[Variables]
    [./disp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_y]
        order = SECOND
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    [./vel_x]
        order = SECOND
        family = LAGRANGE
    []
    [./accel_x]
    []
    [./vel_y]
    []
    [./accel_y]
    []
    [./nodal_area]
        order = SECOND
        family = LAGRANGE
    [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y normal_traction tangent_traction normal_jump tangent_jump'
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
    [./czm_stress_derivative]
        type = StressDerivative2
        boundary = 'Block0_Block1'
        args = 'disp_x disp_y' 
    [../]
    [./czm_mat]
        type = SlipWeakeningFriction2dxx
        boundary = 'Block0_Block1'
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
   # solve_type = 'PJFNK'
    dt = 0.0025
    end_time = 12
    automatic_scaling = true
    [TimeIntegrator]
        type = CentralDifference
      #  type = NewmarkBeta
        #solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 10
[]
