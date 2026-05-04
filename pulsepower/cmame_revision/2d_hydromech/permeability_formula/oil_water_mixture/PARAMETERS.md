# Oil-fluid parameter choice — `oil_water_mixture/`

This case mirrors the water baseline in `permeability_formula/` but swaps the
single-phase pore fluid from water to a generic light crude oil. The model
itself (Heider-2021 normal-strain phase-field hydraulic fracture, fully
saturated single-component PorousFlow) is unchanged — only the three fluid
properties consumed by `SimpleFluidProperties` differ.

## Reference

| Source | Use |
|---|---|
| Batzle, M. & Wang, Z. (1992). "Seismic properties of pore fluids." *Geophysics* 57(11), 1396–1408. [DOI: 10.1190/1.1443207](https://library.seg.org/doi/10.1190/1.1443207) | Closed-form correlations for ρ, K, μ of oil at reservoir T,P. |
| Feng, Y., Haugen, K., Firoozabadi, A. (2021). "Phase-Field Simulation of Hydraulic Fracturing by CO₂, Water and Nitrogen in 2D and Comparison With Laboratory Data." *JGR: Solid Earth* 126, e2021JB022509. [DOI: 10.1029/2021JB022509](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB022509) | Methodological precedent: same single-fluid phase-field model run with three different fluids by changing only ρ, K, μ. |

## Reservoir state assumed

A "generic light crude" at typical reservoir conditions:

| quantity | value |
|---|---|
| API gravity | ~30° (light crude) |
| Reservoir temperature, T | ~80 °C |
| Reservoir pressure, P | ~20 MPa |
| Gas-oil ratio (GOR), R_g | ~0 (dead oil baseline) |

These are within the range Batzle & Wang (1992) explicitly tabulate and
within the parameter envelope for which their correlations were validated.

## Fluid parameter values

These three scalars are defined at the top of `static_solve.i` and
`elasticity_E1d25.i`; every downstream block (`SimpleFluidProperties`,
`ElkPorousFlowDamagedBiotModulus`, kinetic-energy postprocessors) reads them
through `${...}` substitution, so editing the header is the *only* change
required to retarget the fluid.

| parameter | water (baseline) | oil (this case) | source |
|---|---|---|---|
| `fluid_density`      | 1000 kg/m³  | **850 kg/m³**  | Batzle-Wang ρ_oil at 30° API, dead-oil baseline. |
| `fluid_bulk_modulus` | 2.24e9 Pa   | **1.5e9 Pa**   | Batzle-Wang K_oil from V_p,oil(API,T,P); ~30% softer than water. |
| `viscosity`          | 1e-3 Pa·s   | **5e-2 Pa·s**  | Batzle-Wang μ_oil for light crude at reservoir T (50× water). |

The viscosity is the parameter with the widest spread across published light-crude
correlations (1e-2 .. 1e-1 Pa·s depending on T and dissolved-gas content). 5e-2
sits near the middle and is a defensible single-value choice; bump up for a
heavier/colder oil, bump down for a hotter oil with significant dissolved gas.

## What is *not* changed (and why)

* **Solid properties** (E, ν, ρ_s, G_c, ℓ): same as the water case. Oil only
  enters the model through the pore fluid, not the rock matrix.
* **Permeability model**: `permeability_model = normal_strain` (Heider 2021,
  eqs. 46-48), `intrinsic_permeability = 5e-19` m², `perm_exponent = 10`. The
  effective permeability is a function of the rock fabric and crack geometry,
  not the fluid; identical to the water case.
* **Biot coupling**: `K_s = 50e9` Pa for the solid grains; the damaged
  Biot coefficient/modulus are computed from K, K_s and the fluid bulk
  modulus. Only the fluid term changes.
* **Boundary conditions**: same `func_tri_pulse` mechanical pressure load on
  the borehole, same confinement (1 MPa), same dampers — the fracturing
  source is mechanical, not chemical, so swapping fluid does not change the
  driving force.

## Files in this directory

```
oil_water_mixture/
├── PARAMETERS.md          # this file
├── static_solve.i         # quasi-static initialisation; oil fluid props
├── elasticity_E1d25.i     # dynamic phase-field run; reads ./static_solve_out.e
└── fracture_E1d25.i       # AT1 phase-field sub-app (no fluid props; unchanged)
```

## Run order

1. `static_solve.i` — produces `static_solve_out.e` in this directory.
2. `elasticity_E1d25.i` — `SolutionUserObject` reads `./static_solve_out.e`
   for `disp_x`, `disp_y`, `pp` initial conditions; then runs the dynamic
   pulse-power loading with the oil pore fluid.

## Citing this choice in a manuscript

A defensible single sentence:

> The oil pore-fluid properties (ρ = 850 kg/m³, K = 1.5 GPa, μ = 5×10⁻² Pa·s)
> were taken from the Batzle & Wang (1992) correlations for a 30° API light
> crude at reservoir conditions (T ≈ 80 °C, P ≈ 20 MPa, R_g ≈ 0), following
> the parametric-fluid-swap approach used by Feng et al. (2021) for CO₂,
> water, and nitrogen.
