#Test the analytical pressure solution
[Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 2
      nx = 100
      ny = 100
      xmin = -2500
      xmax = 2500
      ymin = -2500
      ymax = 2500
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'

    #pressure analytical solution
    #Reference: Injection-induced seismicity: Poroelastic and earthquake nucleation effects (P. Segall1 and S. Lu2)
    # flux_q = 1e-2 #kg/s
    # density_rho_0 = 1e3 #kg/m^3
    # permeability_k = 3e-16 #m^2
    # viscosity_eta = 0.4e-3 #Pa s
    # biotcoeff_alpha = 0.31 #-
    # undrained_nu_u = 0.3  #-
    # shear_modulus_mu = 20e9 #Pa
    # drained_nu = 0.25 #-
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
        [../]
      [../]
    [../]
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
[]

[AuxVariables]
    [pressure]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxKernels]
    #initial is zero, DO NOT add initial on execute_on
    [./aux_func_pressure]
        type = FunctionAux
        variable = pressure
        function = func_pressure
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
[]

[Functions]
    [func_pressure]
        type = PiecewiseMultilinear
        data_file = analyticalsol.txt
    []
[]

[Executioner]
    type = Transient
    dt = 0.1
    end_time = 10.0
    # num_steps = 10
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 10
[]