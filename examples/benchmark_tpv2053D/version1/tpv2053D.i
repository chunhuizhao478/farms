# Verification of Benchmark Problem TPV205-3D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 150
    ny = 150
    nz = 150
    xmin = -15000
    xmax = 15000
    ymin = -15000
    ymax = 15000
    zmin = -15000
    zmax = 15000
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
    displacements = 'disp_x disp_y disp_z'
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
    area = 200
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
    [./vel_z]
    []
    [./accel_z]
    []
[]
  
[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
      boundary = 'Block0_Block1'
      strain = SMALL
    [../]
[]
  
[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          add_variables = true
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
    [./inertia_z]
      type = InertialForce
      use_displaced_mesh = false
      variable = disp_z
    []
    [./StressDampingx]
      type = StiffPropDamping
      variable = 'disp_x'
      component = '0'
    []
    [./StressDampingy]
      type = StiffPropDamping
      variable = 'disp_y'
      component = '1'
    []
    [./StressDampingz]
      type = StiffPropDamping
      variable = 'disp_z'
      component = '2'
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
        type = SlipWeakeningFriction3d
        boundary = 'Block0_Block1'
    [../]
[]
  
[Executioner]
    type = Transient
    dt = 0.005
    end_time = 3
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
[]
  
[Outputs]
    exodus = true
    interval = 10
[]