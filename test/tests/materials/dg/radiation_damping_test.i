# Radiation Damping Material Test
# Verify radiation damping coefficient calculation for BP2 parameters

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[Materials]
  [rad_damp]
    type = RadiationDampingMaterial
    shear_modulus = 32.04e9   # Pa (BP2)
    density = 2670            # kg/m³ (BP2)
  []
[]

[Executioner]
  type = Steady
[]

[Postprocessors]
  # Expected: cs = sqrt(32.04e9 / 2670) ≈ 3464 m/s
  [shear_wave_speed]
    type = ElementAverageMaterialProperty
    mat_prop = shear_wave_speed
  []
  # Expected: eta = 32.04e9 / (2 * 3464) ≈ 4.625e6 Pa·s/m
  [radiation_damping]
    type = ElementAverageMaterialProperty
    mat_prop = radiation_damping
  []
[]

[Outputs]
  csv = true
[]
