# =============================================================================
# SPE-HFTC acrylic cube, pulse-power fracturing, pure solid, 3D (main app)
# ABSORBING-FACE VARIANT: Lysmer dashpots on the four side faces (xmin xmax ymin ymax) instead of free faces;
# the end faces keep the holder condition u_z = 0. Everything else identical to elasticity.i; the dashpot work
# enters the energy budget as damping_work. This is the 3D twin of ../code_2d/elasticity_absorbing.i.
# Dynamic elasticity (Newmark) + AT1 phase-field fracture in the sub-app fracture.i
# Decisions of 2026-10-21 are marked [D]; see ../PLAN.md and ../EXPLORE_3d_solid.md.
# Mesh: ../mesh/make_spe_cube_mesh.py (named sidesets borehole xmin xmax ymin ymax zmin zmax, block acrylic)
# Local run (Mac):  PATH=../tools/jit_shim:$PATH mpiexec -n 8 farms-opt -i elasticity_absorbing.i mesh_file=../mesh/spe_cube_smoke.msh l=3e-3 Np=1
# Vista (production mesh, 14.6 M tets): split first, see cases_vista/README.md
# =============================================================================

# ---- material: cast acrylic, paper Table 1 (US Plastic Corp. cast acrylic sheet GP) ----
E = 2.80e9
nu = 0.35
density = 1190
ft = 69e6                          # tensile strength (ASTM D638)
l = 6e-4                          # phase-field length: 2 elements of the default mesh
Gc_const = ${fparse 8*l*ft*ft/(3*E)}   # [D] AT1 strength relation inverted: sigma_c = ft  =>  Gc = 8 l ft^2/(3E)
                                   #     = 2721 J/m^2 at l = 0.6 mm (7.8x the literature 350 J/m^2 of PMMA)
K = ${fparse E/3.0/(1.0-2.0*nu)}
G = ${fparse E/2.0/(1.0+nu)}
Cs = ${fparse sqrt(G/density)}                 # dashpot wave speeds
Cp = ${fparse sqrt((K + 4.0/3.0*G)/density)}

# ---- geometry: paper Table 2 ----
Lx = 51.11e-3
Ly = 51.92e-3
Lz = 49.14e-3
rb = 4.15e-3                       # [D] 8.3 mm borehole (Fig. 2 caption; Table 2 says 7.94 mm)
gap = 7.17e-3                      # electrode gap (Table 2)
zc = ${fparse Lz/3}                # [D] gap centred at one third of the height, measured from the bottom face
mesh_file = ../mesh/spe_cube_h0p3.msh

# ---- pulse: paper Table 2 + empirical peak-pressure law coded in ElkPulseLoadExperiment ----
EM = 0.144                         # stored energy per pulse, kJ (60 kV, 80 nF)
eta = 0.3                          # [D] fraction of EM reaching the wall: free parameter, swept 0.1-1
base_factor = 8000                 # [D] peak-pressure law constants kept from the granite cases (bar, mm, kJ)
kappa = 0.35
m_exp = 0.25
t_r = 3e-7                         # [D] pressure rise time = current rise time (230-350 ns)
# discharge duration to 10 %: paper 3-12 us. Shape rates (peak at t_r, 10 % at t_dur):
#   t_dur =  3 us: shape_alpha = 8.9489e5  shape_beta = 8.3321e6
#   t_dur =  6 us: shape_alpha = 4.1032e5  shape_beta = 1.1529e7   (default)
#   t_dur = 12 us: shape_alpha = 1.9798e5  shape_beta = 1.4514e7
shape_alpha = 4.1032e5
shape_beta = 1.1529e7
Tw = 6e-5                          # window per pulse: pulse + ring-down (quench)
Np = 20                            # [D] pulses in one continuous train, no reset between pulses; count not binding
                                   #     (the time-point lists below cover Np <= 40)

# ---- damping between pulses (decision pending; recommended quench, see PLAN.md Appendix B) ----
t_free = 2.5e-5                    # undamped part of each window: pulse, primary crack growth, first reflection
zeta_pulse = 0                     # Rayleigh stiffness damping (s) during t_free
zeta_quench = 1e-5                 # Rayleigh stiffness damping (s) after t_free: brings the cube to rest
                                   # constant damping: zeta_pulse = zeta_quench ; none: both 0
hht_alpha = 0

# ---- time stepping ----
dt0 = 1e-8                         # step at each pulse onset
dt_rise = 3e-8                     # cap during the first t_rise_hold of each window (10 steps on the 0.3 us rise)
t_rise_hold = 6e-7
dt_max = 1e-7
dt_quiet = ${dt_max}                # cap after t_free of each window (ring-down under zeta_quench); must be >= dt_max;
                                   # the 2D runs so far used dt_quiet = dt_max = 1e-7

# ---- metrics ----
d_th = 0.5                         # damage threshold for the fracture extent (d >= 0.9 also reported)
r_core = 15e-3                     # radius of the "core" region: extent/area reported inside it exclude face-spall and corner damage

# note: the AD parsed materials (degradation, crack geometric function, psi) require the FParser JIT;
#       on this Mac the JIT linker fails (see memory note toolchain-macos-tahoe-link-failure) unless
#       FPARSER_JIT_* flags are set, see ../mesh/README.md
newmark_beta = 0.25
newmark_gamma = 0.5

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc_const=${Gc_const};l=${l};mesh_file=${mesh_file}'
    execute_on = 'INITIAL TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = 'd dissipated_energy_density Gc_var'
    source_variable = 'd dissipated_energy_density Gc_var'
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active mesh_size'
  []
  [pp_transfer_dissipated_energy_total]
    type = MultiAppPostprocessorTransfer
    from_multi_app = 'fracture'
    from_postprocessor = 'dissipated_energy_dynamic'
    to_postprocessor = 'dissipated_energy_dynamic'
    reduction_type = 'sum'
  []
  [pp_transfer_dissipated_energy_first_step]
    type = MultiAppPostprocessorTransfer
    from_multi_app = 'fracture'
    from_postprocessor = 'dissipated_energy_first_step'
    to_postprocessor = 'dissipated_energy_first_step'
    reduction_type = 'sum'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [msh]
    type = FileMeshGenerator
    file = ${mesh_file}
  []
  [corner]
    type = ExtraNodesetGenerator
    input = msh
    coord = '${fparse -Lx/2} ${fparse -Ly/2} 0'
    new_boundary = corner_ptr
    use_closest_node = true
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = FIRST
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
  []
  [disp_z]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [eng_err]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  [vel_z]
    family = LAGRANGE
    order = FIRST
  []
  [accel_z]
    family = LAGRANGE
    order = FIRST
  []
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  [dissipated_energy_density]
    family = MONOMIAL
    order = CONSTANT
  []
  # reaction forces of the borehole pressure (save_in), used for the external work
  [fx]
  []
  [fy]
  []
  [fz]
  []
  # dashpot reaction forces (save_in), for the damping work
  [fdampx]
  []
  [fdampy]
  []
  [fdampz]
  []
  # fracture metrics
  [r_damaged]
    family = LAGRANGE
    order = FIRST
  []
  [r_damaged_090]
    family = LAGRANGE
    order = FIRST
  []
  [r_damaged_core]
    family = LAGRANGE
    order = FIRST
  []
  [is_damaged_core]
    family = MONOMIAL
    order = CONSTANT
  []
  [is_damaged]
    family = MONOMIAL
    order = CONSTANT
  []
  [solid_kinetic_energy]
    order = CONSTANT
    family = MONOMIAL
  []
  # particle velocity magnitude (documented in the Exodus output and as v_max in the CSV)
  [vel_mag]
    family = LAGRANGE
    order = FIRST
  []
  # toughness field of the sub-app (constant, Weibull or Gaussian), copied here for output
  [Gc_var]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [error_measure]
    type = ErrorPsiMeasure
    variable = eng_err
  []
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = ${newmark_beta}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = ${newmark_gamma}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_pulse_load_aux]
    type = FunctionAux
    variable = pulse_load_aux
    function = func_pulse
    execute_on = timestep_end
  []
  [mesh_size_aux]
    type = MeshSize
    variable = mesh_size
    execute_on = 'INITIAL'
  []
  # radial distance (from the borehole axis) of every node with d >= threshold, 0 elsewhere
  [r_damaged_aux]
    type = ParsedAux
    variable = r_damaged
    coupled_variables = 'd'
    use_xyzt = true
    expression = 'if(d >= ${d_th}, sqrt(x*x + y*y), 0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [r_damaged_090_aux]
    type = ParsedAux
    variable = r_damaged_090
    coupled_variables = 'd'
    use_xyzt = true
    expression = 'if(d >= 0.9, sqrt(x*x + y*y), 0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [is_damaged_aux]
    type = ParsedAux
    variable = is_damaged
    coupled_variables = 'd'
    expression = 'if(d >= ${d_th}, 1, 0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # same two measures restricted to r <= r_core around the borehole (radial cracks only)
  [r_damaged_core_aux]
    type = ParsedAux
    variable = r_damaged_core
    coupled_variables = 'd'
    use_xyzt = true
    expression = 'if(d >= ${d_th} & sqrt(x*x + y*y) <= ${r_core}, sqrt(x*x + y*y), 0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [is_damaged_core_aux]
    type = ParsedAux
    variable = is_damaged_core
    coupled_variables = 'd'
    use_xyzt = true
    expression = 'if(d >= ${d_th} & sqrt(x*x + y*y) <= ${r_core}, 1, 0)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_mag_aux]
    type = ParsedAux
    variable = vel_mag
    coupled_variables = 'vel_x vel_y vel_z'
    expression = 'sqrt(vel_x*vel_x + vel_y*vel_y + vel_z*vel_z)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [solid_kinetic_energy]
    type = ADKineticEnergyAux
    variable = solid_kinetic_energy
    newmark_velocity_x = vel_x
    newmark_velocity_y = vel_y
    newmark_velocity_z = vel_z
    density = density
  []
[]

[Kernels]
  [solid_x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    zeta = zeta                      # material property: time-switched Rayleigh damping (zeta_fn)
    displacements = 'disp_x disp_y disp_z'
  []
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_x
    acceleration = accel_x
    density = ${density}
    alpha = ${hht_alpha}
  []
  [solid_y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    zeta = zeta                      # material property: time-switched Rayleigh damping (zeta_fn)
    displacements = 'disp_x disp_y disp_z'
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_y
    acceleration = accel_y
    density = ${density}
    alpha = ${hht_alpha}
  []
  [solid_z]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = false
    alpha = ${hht_alpha}
    zeta = zeta                      # material property: time-switched Rayleigh damping (zeta_fn)
    displacements = 'disp_x disp_y disp_z'
  []
  [inertia_z]
    type = ADInertialForce
    variable = disp_z
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    velocity = vel_z
    acceleration = accel_z
    density = ${density}
    alpha = ${hht_alpha}
  []
[]

[Functions]
  # borehole pressure p(t, x): double-exponential pulse repeated every Tw, peak from the empirical law
  [func_pulse]
    type = ElkPulseLoadExperiment
    shape_param_alpha = ${shape_alpha}
    shape_param_beta = ${shape_beta}
    rise_time = ${t_r}
    single_pulse_duration = ${Tw}
    EM = ${EM}
    gap = ${gap}
    convert_efficiency = ${eta}
    fitting_param_alpha = ${kappa}
    fitting_param_exponent = ${m_exp}
    base_factor = ${base_factor}
    discharge_center = '0 0 ${zc}'      # gap centred at zc on the borehole axis
    number_of_pulses = ${Np}
  []
  # pulse onsets k*Tw: the time stepper lands exactly on them and restarts at dt0
  [pulse_sync]
    type = PiecewiseLinear
    x = '${fparse 0*Tw} ${fparse 1*Tw} ${fparse 2*Tw} ${fparse 3*Tw} ${fparse 4*Tw} ${fparse 5*Tw} ${fparse 6*Tw} ${fparse 7*Tw} ${fparse 8*Tw} ${fparse 9*Tw} ${fparse 10*Tw} ${fparse 11*Tw} ${fparse 12*Tw} ${fparse 13*Tw} ${fparse 14*Tw} ${fparse 15*Tw} ${fparse 16*Tw} ${fparse 17*Tw} ${fparse 18*Tw} ${fparse 19*Tw} ${fparse 20*Tw} ${fparse 21*Tw} ${fparse 22*Tw} ${fparse 23*Tw} ${fparse 24*Tw} ${fparse 25*Tw} ${fparse 26*Tw} ${fparse 27*Tw} ${fparse 28*Tw} ${fparse 29*Tw} ${fparse 30*Tw} ${fparse 31*Tw} ${fparse 32*Tw} ${fparse 33*Tw} ${fparse 34*Tw} ${fparse 35*Tw} ${fparse 36*Tw} ${fparse 37*Tw} ${fparse 38*Tw} ${fparse 39*Tw} ${fparse 40*Tw}'
    y = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40'
  []
  # time-step cap: dt_rise for the first t_rise_hold of each window, dt_max until t_free (pulse and crack growth),
  # dt_quiet in the damped ring-down after t_free (dt_quiet = dt_max reproduces the schedule used so far)
  [dt_limit]
    type = PiecewiseConstant
    direction = LEFT
    x = '${fparse 0*Tw} ${fparse 0*Tw + t_rise_hold} ${fparse 0*Tw + t_free} ${fparse 1*Tw} ${fparse 1*Tw + t_rise_hold} ${fparse 1*Tw + t_free} ${fparse 2*Tw} ${fparse 2*Tw + t_rise_hold} ${fparse 2*Tw + t_free} ${fparse 3*Tw} ${fparse 3*Tw + t_rise_hold} ${fparse 3*Tw + t_free} ${fparse 4*Tw} ${fparse 4*Tw + t_rise_hold} ${fparse 4*Tw + t_free} ${fparse 5*Tw} ${fparse 5*Tw + t_rise_hold} ${fparse 5*Tw + t_free} ${fparse 6*Tw} ${fparse 6*Tw + t_rise_hold} ${fparse 6*Tw + t_free} ${fparse 7*Tw} ${fparse 7*Tw + t_rise_hold} ${fparse 7*Tw + t_free} ${fparse 8*Tw} ${fparse 8*Tw + t_rise_hold} ${fparse 8*Tw + t_free} ${fparse 9*Tw} ${fparse 9*Tw + t_rise_hold} ${fparse 9*Tw + t_free} ${fparse 10*Tw} ${fparse 10*Tw + t_rise_hold} ${fparse 10*Tw + t_free} ${fparse 11*Tw} ${fparse 11*Tw + t_rise_hold} ${fparse 11*Tw + t_free} ${fparse 12*Tw} ${fparse 12*Tw + t_rise_hold} ${fparse 12*Tw + t_free} ${fparse 13*Tw} ${fparse 13*Tw + t_rise_hold} ${fparse 13*Tw + t_free} ${fparse 14*Tw} ${fparse 14*Tw + t_rise_hold} ${fparse 14*Tw + t_free} ${fparse 15*Tw} ${fparse 15*Tw + t_rise_hold} ${fparse 15*Tw + t_free} ${fparse 16*Tw} ${fparse 16*Tw + t_rise_hold} ${fparse 16*Tw + t_free} ${fparse 17*Tw} ${fparse 17*Tw + t_rise_hold} ${fparse 17*Tw + t_free} ${fparse 18*Tw} ${fparse 18*Tw + t_rise_hold} ${fparse 18*Tw + t_free} ${fparse 19*Tw} ${fparse 19*Tw + t_rise_hold} ${fparse 19*Tw + t_free} ${fparse 20*Tw} ${fparse 20*Tw + t_rise_hold} ${fparse 20*Tw + t_free} ${fparse 21*Tw} ${fparse 21*Tw + t_rise_hold} ${fparse 21*Tw + t_free} ${fparse 22*Tw} ${fparse 22*Tw + t_rise_hold} ${fparse 22*Tw + t_free} ${fparse 23*Tw} ${fparse 23*Tw + t_rise_hold} ${fparse 23*Tw + t_free} ${fparse 24*Tw} ${fparse 24*Tw + t_rise_hold} ${fparse 24*Tw + t_free} ${fparse 25*Tw} ${fparse 25*Tw + t_rise_hold} ${fparse 25*Tw + t_free} ${fparse 26*Tw} ${fparse 26*Tw + t_rise_hold} ${fparse 26*Tw + t_free} ${fparse 27*Tw} ${fparse 27*Tw + t_rise_hold} ${fparse 27*Tw + t_free} ${fparse 28*Tw} ${fparse 28*Tw + t_rise_hold} ${fparse 28*Tw + t_free} ${fparse 29*Tw} ${fparse 29*Tw + t_rise_hold} ${fparse 29*Tw + t_free} ${fparse 30*Tw} ${fparse 30*Tw + t_rise_hold} ${fparse 30*Tw + t_free} ${fparse 31*Tw} ${fparse 31*Tw + t_rise_hold} ${fparse 31*Tw + t_free} ${fparse 32*Tw} ${fparse 32*Tw + t_rise_hold} ${fparse 32*Tw + t_free} ${fparse 33*Tw} ${fparse 33*Tw + t_rise_hold} ${fparse 33*Tw + t_free} ${fparse 34*Tw} ${fparse 34*Tw + t_rise_hold} ${fparse 34*Tw + t_free} ${fparse 35*Tw} ${fparse 35*Tw + t_rise_hold} ${fparse 35*Tw + t_free} ${fparse 36*Tw} ${fparse 36*Tw + t_rise_hold} ${fparse 36*Tw + t_free} ${fparse 37*Tw} ${fparse 37*Tw + t_rise_hold} ${fparse 37*Tw + t_free} ${fparse 38*Tw} ${fparse 38*Tw + t_rise_hold} ${fparse 38*Tw + t_free} ${fparse 39*Tw} ${fparse 39*Tw + t_rise_hold} ${fparse 39*Tw + t_free}'
    y = '${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet} ${dt_rise} ${dt_max} ${dt_quiet}'
  []
  # Rayleigh damping schedule: zeta_pulse for the first t_free of each window, zeta_quench afterwards
  [zeta_fn]
    type = PiecewiseConstant
    direction = LEFT
    x = '${fparse 0*Tw} ${fparse 0*Tw + t_free} ${fparse 1*Tw} ${fparse 1*Tw + t_free} ${fparse 2*Tw} ${fparse 2*Tw + t_free} ${fparse 3*Tw} ${fparse 3*Tw + t_free} ${fparse 4*Tw} ${fparse 4*Tw + t_free} ${fparse 5*Tw} ${fparse 5*Tw + t_free} ${fparse 6*Tw} ${fparse 6*Tw + t_free} ${fparse 7*Tw} ${fparse 7*Tw + t_free} ${fparse 8*Tw} ${fparse 8*Tw + t_free} ${fparse 9*Tw} ${fparse 9*Tw + t_free} ${fparse 10*Tw} ${fparse 10*Tw + t_free} ${fparse 11*Tw} ${fparse 11*Tw + t_free} ${fparse 12*Tw} ${fparse 12*Tw + t_free} ${fparse 13*Tw} ${fparse 13*Tw + t_free} ${fparse 14*Tw} ${fparse 14*Tw + t_free} ${fparse 15*Tw} ${fparse 15*Tw + t_free} ${fparse 16*Tw} ${fparse 16*Tw + t_free} ${fparse 17*Tw} ${fparse 17*Tw + t_free} ${fparse 18*Tw} ${fparse 18*Tw + t_free} ${fparse 19*Tw} ${fparse 19*Tw + t_free} ${fparse 20*Tw} ${fparse 20*Tw + t_free} ${fparse 21*Tw} ${fparse 21*Tw + t_free} ${fparse 22*Tw} ${fparse 22*Tw + t_free} ${fparse 23*Tw} ${fparse 23*Tw + t_free} ${fparse 24*Tw} ${fparse 24*Tw + t_free} ${fparse 25*Tw} ${fparse 25*Tw + t_free} ${fparse 26*Tw} ${fparse 26*Tw + t_free} ${fparse 27*Tw} ${fparse 27*Tw + t_free} ${fparse 28*Tw} ${fparse 28*Tw + t_free} ${fparse 29*Tw} ${fparse 29*Tw + t_free} ${fparse 30*Tw} ${fparse 30*Tw + t_free} ${fparse 31*Tw} ${fparse 31*Tw + t_free} ${fparse 32*Tw} ${fparse 32*Tw + t_free} ${fparse 33*Tw} ${fparse 33*Tw + t_free} ${fparse 34*Tw} ${fparse 34*Tw + t_free} ${fparse 35*Tw} ${fparse 35*Tw + t_free} ${fparse 36*Tw} ${fparse 36*Tw + t_free} ${fparse 37*Tw} ${fparse 37*Tw + t_free} ${fparse 38*Tw} ${fparse 38*Tw + t_free} ${fparse 39*Tw} ${fparse 39*Tw + t_free} ${fparse 40*Tw} ${fparse 40*Tw + t_free}'
    y = '${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench} ${zeta_pulse} ${zeta_quench}'
  []
[]

[BCs]
  [./Pressure]
    [pressure_borehole]
      boundary = borehole
      function = func_pulse
      displacements = 'disp_x disp_y disp_z'
      use_displaced_mesh = false
      save_in_disp_x = fx
      save_in_disp_y = fy
      save_in_disp_z = fz
    []
  []
  # absorbing side faces: Lysmer dashpots (rho c_p v_n in the normal, rho c_s v_t in the tangential directions)
  [damp_x]
    type = FarmsNonReflectDashpotBC
    variable = disp_x
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 0
    boundary = 'xmin xmax ymin ymax'
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampx
  []
  [damp_y]
    type = FarmsNonReflectDashpotBC
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 1
    boundary = 'xmin xmax ymin ymax'
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampy
  []
  [damp_z]
    type = FarmsNonReflectDashpotBC
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    component = 2
    boundary = 'xmin xmax ymin ymax'
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    shear_wave_speed = ${Cs}
    p_wave_speed = ${Cp}
    density = ${density}
    save_in = fdampz
  []
  # [D] holder plates (paper Fig. 2a): no vertical motion of the end faces (the side faces carry the dashpots)
  [fix_z_ends]
    type = DirichletBC
    variable = disp_z
    boundary = 'zmin zmax'
    value = 0
  []
  # one corner node pinned in-plane (removes in-plane rigid modes; harmless for the dynamics)
  [fix_corner_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr
    value = 0
  []
  [fix_corner_y]
    type = DirichletBC
    variable = disp_y
    boundary = corner_ptr
    value = 0
  []
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G density'
    prop_values = '${K} ${G} ${density}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [zeta_mat]
    type = GenericFunctionMaterial
    prop_names = 'zeta'
    prop_values = 'zeta_fn'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -snes_type -ksp_gmres_restart'
  petsc_options_value = 'gmres hypre boomeramg 0.7 4 5 0.3 vinewtonrsls 100'

  line_search = 'basic'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 30

  end_time = ${fparse Np*Tw}

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8

  [TimeStepper]
    type = FarmsIterationAdaptiveDT
    dt = ${dt0}
    iteration_window = 0
    cutback_factor_at_failure = 0.5
    optimal_iterations = 20
    growth_factor = 1.25
    max_time_step_bound = ${dt_quiet}
    timestep_limiting_function = pulse_sync
    force_step_every_function_point = true
    post_function_sync_dt = ${dt0}
    timestep_limiting_postprocessor = dt_limit_pp
  []
  [TimeIntegrator]
    type = NewmarkBeta
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    time_step_interval = 100       # 14.6 M tets: about 6 frames per pulse, nodal fields only (19 MB per field and frame)
    execute_on = 'INITIAL TIMESTEP_END FINAL'
    show = 'd vel_mag'
  []
  [checkpoint]
    type = Checkpoint
    time_step_interval = 100
    num_files = 2
  []
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
    time_step_interval = 1
    show = 'full_energy full_input_energy external_work damping_work solid_elastic_energy_total solid_kinetic_energy_total dissipated_energy_total fracture_extent fracture_extent_090 fracture_extent_core r_dmg_max r_dmg_max_090 damaged_volume damaged_volume_core crack_area p_wall_max v_max borehole_area'
  []
  [perf]
    # wall-time breakdown (solve, sub-app, transfers, aux, postprocessors, output) printed to the log every 500 steps,
    # so that a segment killed by the wall limit still leaves its cost profile; the timing runs override the interval
    type = PerfGraphOutput
    execute_on = 'TIMESTEP_END FINAL'
    time_step_interval = 500
    level = 2
  []
[]

# ============================ energy bookkeeping =============================
[Postprocessors]
  # crack energy, integrated in the sub-app and received here; the first-step value is subtracted
  [dissipated_energy_dynamic]
    type = Receiver
  []
  [dissipated_energy_first_step]
    type = Receiver
  []
  [dissipated_energy_total]
    type = ParsedPostprocessor
    expression = 'dissipated_energy_dynamic - dissipated_energy_first_step'
    pp_names = 'dissipated_energy_dynamic dissipated_energy_first_step'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # work of the borehole pressure (nodal, on the pressure sideset); the clamped ends do no work
  [external_work]
    type = FarmsExternalWork
    boundary = borehole
    forces = 'fx fy fz'
  []
  # work absorbed by the dashpots on the side faces
  [damping_work]
    type = FarmsExternalWork
    boundary = 'xmin xmax ymin ymax'
    forces = 'fdampx fdampy fdampz'
  []
  [full_input_energy]
    type = ParsedPostprocessor
    expression = '-1 * external_work - damping_work'
    pp_names = 'external_work damping_work'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [solid_kinetic_energy_total]
    type = ElementIntegralVariablePostprocessor
    variable = solid_kinetic_energy
  []
  [solid_elastic_energy_dynamic]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie
  []
  [solid_elastic_energy_total]
    type = ParsedPostprocessor
    expression = 'solid_elastic_energy_dynamic'
    pp_names = 'solid_elastic_energy_dynamic'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # budget check: full_energy = full_input_energy when zeta = 0 (quench energy is the gap otherwise)
  [full_energy]
    type = ParsedPostprocessor
    expression = 'solid_kinetic_energy_total + solid_elastic_energy_total + dissipated_energy_total'
    pp_names = 'solid_kinetic_energy_total solid_elastic_energy_total dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # time-step cap schedule (read by the time stepper)
  [dt_limit_pp]
    type = FunctionValuePostprocessor
    function = dt_limit
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# ============================ fracture metrics ===============================
[Postprocessors]
  # longest fracture from the borehole wall in top view: max radial distance of nodes with d >= d_th, minus rb
  [r_dmg_max]
    type = NodalExtremeValue
    variable = r_damaged
    value_type = max
  []
  [fracture_extent]
    type = ParsedPostprocessor
    expression = 'max(r_dmg_max - ${rb}, 0)'
    pp_names = 'r_dmg_max'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [r_dmg_max_090]
    type = NodalExtremeValue
    variable = r_damaged_090
    value_type = max
  []
  [fracture_extent_090]
    type = ParsedPostprocessor
    expression = 'max(r_dmg_max_090 - ${rb}, 0)'
    pp_names = 'r_dmg_max_090'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [r_dmg_max_core]
    type = NodalExtremeValue
    variable = r_damaged_core
    value_type = max
  []
  [fracture_extent_core]
    type = ParsedPostprocessor
    expression = 'max(r_dmg_max_core - ${rb}, 0)'
    pp_names = 'r_dmg_max_core'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damaged_volume_core]
    type = ElementIntegralVariablePostprocessor
    variable = is_damaged_core
  []
  # damaged volume (d >= d_th), m^3
  [damaged_volume]
    type = ElementIntegralVariablePostprocessor
    variable = is_damaged
  []
  # crack area estimate = crack energy / Gc, m^2 (divide by 7.8 to compare with a 350 J/m^2 material)
  [crack_area]
    type = ParsedPostprocessor
    expression = 'dissipated_energy_total / ${Gc_const}'
    pp_names = 'dissipated_energy_total'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # peak of the applied pressure (the function's maximum over the solid sits on the borehole wall at the gap)
  [p_wall_max]
    type = ElementExtremeValue
    variable = pulse_load_aux
    value_type = max
  []
  # sideset check: 2*pi*rb*Lz = 1.2813e-3 m^2
  # maximum particle velocity over the domain
  [v_max]
    type = NodalExtremeValue
    variable = vel_mag
    value_type = max
  []
  [borehole_area]
    type = AreaPostprocessor
    boundary = borehole
    execute_on = INITIAL
  []
[]
