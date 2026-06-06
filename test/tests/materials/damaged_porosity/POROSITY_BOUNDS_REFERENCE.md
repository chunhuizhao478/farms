# Reference & justification — porosity bounds in `ElkPorousFlowDamagedPorosity`

This document records the physical and numerical justification for the porosity
clamp bounds used by the damage-dependent porosity material
`ElkPorousFlowDamagedPorosity` (and its AD twin `ElkADPorousFlowDamagedPorosity`),
in particular the production values

```
porosity_lower_bound = 0.008
porosity_upper_bound = 0.999
```

It follows the precedent of
`test/tests/materials/permeability_normal_strain/EXPECTED_VALUES.md`: a written
record so a future reviewer can confirm the bounds encode physics, not an
undocumented numerical fudge.

**Bottom line up front.** `0.999` is **not** a measured granite porosity — no
saturated granite, however fractured, comes close (measured maxima are ~6.5 %,
§3). `0.999` is the **numerical regularization of the model's `phi -> 1` void
limit**: a fully phase-field-damaged quadrature point represents an *open,
fluid-filled crack* (a void with a vanishing solid skeleton), and the cap
retains a `1 - 0.999 = 1e-3` solid fraction so the downstream poromechanical
coefficients stay finite (§4). The references below justify (i) the **lower**
bound / `initial_porosity` from intact-granite porosity data, and (ii) the
**upper** bound as the regularized void limit of the phase-field framework the
code already follows. **Note:** the default/production bound is `0.999`; the
separate `porosity_bounded` sensitivity case caps at `0.065` instead — the
measured-max fractured-granite porosity (§3) — for a physically-bounded
granular-skeleton run (§4a).

---

## 1. The model and its limits

Per quadrature point, with damage / phase field `d in [0,1]`
(`src/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.C:41-47`):

```
g(d)    = (1 - d)^2                              # AT1/AT2 degradation
phi_raw = phi_0 + (1 - phi_0) * (1 - g(d))       # = phi_0 + (1 - phi_0)*(1 - (1-d)^2)
phi     = clamp(phi_raw, lower_bound, upper_bound)
```

Limiting behaviour:

| state            | `d` | `g` | `phi_raw`                       |
|------------------|-----|-----|---------------------------------|
| intact           | 0   | 1   | `phi_0`  (= 0.008 in production) |
| fully damaged    | 1   | 0   | `phi_0 + (1 - phi_0) = 1.0`      |

So the **upper bound is a regularized stand-in for the `d -> 1` limit
`phi = 1`**, and the **lower bound equals the intact porosity `phi_0`** (in
production `phi_0 = ${porosity} = 0.008`, so the clamp is a no-op floor at the
intact value).

---

## 2. Lower bound `0.008` — intact granite porosity

Intact, unaltered granite has very low porosity, dominated by microcracks:

- **Schild, Siegesmund, Vollbrecht & Mazurek (2001)** report, for granite
  matrix porosity, **in-situ (impregnated) porosity 0.55–0.59 %**, **laboratory
  (non-impregnated) porosity 1.0–1.17 %** (lab values inflated ×2–2.5 by stress
  release / preparation), and cite earlier work giving **0.8–1.2 vol %** for
  unaltered granite.
  *Geophys. J. Int.* **146**(1), 111–125. DOI
  [10.1046/j.0956-540x.2001.01427.x](https://doi.org/10.1046/j.0956-540x.2001.01427.x).
  **[VERIFIED — fetched]**

- **Westerly granite** (the canonical low-porosity reference granite; 30 % qtz /
  30 % oligoclase / 30 % microcline / 10 % biotite) is widely reported with
  porosity **≈ 1.4 %**. **[VERIFIED via cross-source search synthesis; corroborating, not a single-source fetch]**

- General petrophysics: unfractured crystalline rocks / granite typically
  **< 2 %** (often **< 1 %** in situ). **[VERIFIED — search synthesis of standard references]**

**Conclusion.** `phi_0 = 0.008` (0.8 %) sits at the low end of the documented
unaltered-granite range (0.8–1.2 vol %, Schild et al. 2001) and is consistent
with the in-situ ~0.5 % microcrack porosity. Setting `porosity_lower_bound =
0.008 = phi_0` keeps the *undamaged* porosity at its physical intact value and
prevents the clamp from ever pushing porosity below the rock's true matrix
porosity.

---

## 3. What `0.999` is **not** — maximum *measured* saturated granite porosity

Real granite porosity rises with fracturing and alteration, but only modestly
compared with `0.999`:

- **Staněk & Géraud (2019)** measured fractured/altered Lipnice granite
  (mercury-intrusion porosimetry, 21 specimens, 6 alteration facies): fresh
  granite as low as **0.3 %**, rising to a **maximum of 6.5 %** in the most
  fractured/altered (cavity-bearing) specimen; permeability spans **5 orders of
  magnitude** across the suite.
  *Solid Earth* **10**, 251–274. DOI
  [10.5194/se-10-251-2019](https://doi.org/10.5194/se-10-251-2019).
  **[VERIFIED — fetched]**

- Literature compilations note fractured granite reaching ~10 %, and
  cataclasis / mineral-dissolution secondary porosity up to ~20 % in extreme
  cases. **[search synthesis; order-of-magnitude only — UNVERIFIED single source]**

**Conclusion.** The maximum porosity granite attains in saturated, fully
fractured/altered conditions is **~6.5 % (up to ~20 % in extreme secondary
porosity)** — three orders of magnitude below `0.999`. Therefore `0.999`
**cannot** be justified as "the maximum porosity granite can reach"; it must be
justified as a model void-limit (§4), and this section is the honest record of
that distinction.

---

## 4. Upper bound `0.999` — regularized `phi -> 1` void limit

In a phase-field fracture model the fully damaged region (`d = 1`) is not
"granite with 99.9 % pore space"; it is an **open crack treated as a
fluid-filled void**, i.e. the solid skeleton has been fully degraded and the
pore space locally occupies the element. The interpolation `phi(d)` driving
porosity to `1` at full damage encodes exactly this: at `d = 1` the material
point is all fluid, no skeleton.

- The repository follows the **phase-field-in-porous-media** framework of
  **Heider (2021)**, *A review on phase-field modeling of hydraulic fracturing*,
  *Eng. Fract. Mech.* **253**, 107881. DOI
  [10.1016/j.engfracmech.2021.107881](https://doi.org/10.1016/j.engfracmech.2021.107881).
  This is the same "Heider 2021" cited by the companion permeability material
  (`ElkPorousFlowPermeabilityDamaged`; see
  `permeability_normal_strain/EXPECTED_VALUES.md`), where the fully cracked zone
  is treated as a fluid conduit (parallel-plate / cubic-law permeability
  `k_w = w_h^2/12`). The void-limit treatment of porosity is the mass-balance
  counterpart of that permeability treatment. **[citation VERIFIED — fetched/searched; the specific void-limit reading is the modeling rationale, consistent with this framework]**

- Supporting (Theory of Porous Media phase-field for saturated media):
  **Heider & Markert (2017)**, *A phase-field modeling approach of hydraulic
  fracture in saturated porous media*, Mech. Research Comm. **[cited via the Heider 2021 review — not independently fetched; UNVERIFIED page-level detail]**

### Why cap strictly below 1 (numerical regularization)

Driving `phi` exactly to `1` removes the solid skeleton entirely and makes the
poromechanical coefficients singular or unphysical. Concrete downstream
consumers in this codebase:

- **Damaged Biot modulus**
  (`src/materials/porousflowmatprops/ElkPorousFlowDamagedBiotModulus.C:115-122`):
  ```
  denom = phi/K_fluid + (alpha - phi)/K_grain
  M     = 1 / max(denom, denominator_floor)
  ```
  As `phi -> 1` with Biot coefficient `alpha < 1`, the term `(alpha - phi)/K_grain`
  goes **negative**; the alternate Biot-modulus form noted in the same file,
  `(1 - alpha)(alpha - phi)/K`, likewise becomes negative — an unphysical
  (negative) storage modulus. The code already guards with a `denominator_floor`,
  but keeping `phi <= 0.999` keeps the porosity in the regime where the
  fluid/solid storage partition remains physically meaningful rather than
  relying on the floor alone.

- **Solid-fraction-weighted mass / inertia terms**
  (e.g. `ElkADThreeFieldHistoryEnergyEnhanced.C`), which carry `(1 - phi)` solid
  contributions and `1/phi^2` factors. The retained `1 - phi >= 1e-3` keeps the
  solid mass/stiffness contribution non-vanishing and the factors finite.

The retained solid fraction `1 - 0.999 = 1e-3` is small enough that the
fully-damaged point is mechanically negligible (a crack) yet large enough to
keep every coefficient finite and signed correctly. An upper clamp strictly
below 1 is standard practice for porosity / void interpolation in phase-field
poromechanics for exactly this reason.

**Conclusion.** `porosity_upper_bound = 0.999` is the numerical regularization
of the physically-correct `phi -> 1` void limit of a fully cracked,
fluid-saturated material. It is justified by the phase-field framework the code
already adopts (Heider 2021), not by a measured granite porosity, and the
specific value preserves a `1e-3` solid fraction that keeps the damaged Biot
modulus and solid-fraction-weighted terms finite and physical.

### 4a. The `porosity_bounded` sensitivity case (`0.065`)

The production default is `0.999` (void-limit regularization, above). A separate
sensitivity case —
`pulsepower/cmame_revision/2d_hydromech/permeability_formula/porosity_bounded/elasticity_E1d25.i`,
copied from `permeability_formula/undrained/em_0p005` — instead sets
`porosity_upper_bound = 0.065`, the **maximum measured fractured/altered granite
porosity** (§3, Staněk & Géraud 2019). This deliberately models the fully-damaged
zone as a *still-granular solid skeleton* (~93.5 % solid fraction) rather than an
open void. It accepts the storage-vs-conduit inconsistency flagged above — the
permeability stays an aperture-based cubic-law conduit while the porosity is
capped low — as the explicit purpose of the test: compare the void-limit
(`0.999`) and physically-bounded (`0.065`) runs. Thus `0.065` is the *physically
justified* cap when the damaged material is intended to remain a granular solid,
while `0.999` remains the default when it is intended as an open fluid-filled
crack.

---

## 5. Where the bound is used (codebase audit)

Audit run on branch `porosity-bounded` (from `cdms`):

- **Defaults** (both materials): `porosity_lower_bound = 0.0`,
  `porosity_upper_bound = 0.999`
  (`ElkPorousFlowDamagedPorosity.C:17-20`, `ElkADPorousFlowDamagedPorosity.C:26-29`).
- **Production inputs** (`cmame_revision/2d_hydromech`, `3d_hydromech`,
  `pf_threefield`, hydraulic-fracturing cases): **89** `.i` files set
  `porosity_upper_bound = 0.999`; **83** set `porosity_lower_bound = 0.008`
  (4 use `${porosity}`, which equals `0.008`).
- **Sensitivity case** (this work): `cmame_revision/2d_hydromech/permeability_formula/porosity_bounded/`
  sets `porosity_upper_bound = 0.065` (the physically-bounded cap, §4a), copied
  from `permeability_formula/undrained/em_0p005`. The "89 files set 0.999" count
  above predates this case and excludes it.
- **Non-production outliers** (not part of the production suite, recorded for
  completeness):
  - `pf_code2d_porousflow/parametric_study/archive/confinement/case_1mpa/test1/elasticity.i`
    uses `0.9`, and `.../test2_porosity_0d5/elasticity.i` uses `0.5` — *archived
    sensitivity tests*.
  - `unit_tests/phase2_kernels/test_{ad_poromechanics_coupling_damaged,
    ad_mass_conservation_damaged, ad_dynamic_darcy_flow_phasefield}.i` use
    `porosity_lower_bound = 0.1` — *kernel test fixtures*, deliberately exercising
    a non-default floor.
- **Unit test for this material**:
  `pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i` sets
  `0.008 / 0.999` and samples `d ≈ 0.167, 0.5, 0.833` (max `phi = 0.972`), so the
  upper clamp itself is never exercised by its postprocessors; the clamp at
  `0.999` only engages for `d > 0.968`.

---

## 6. References (verification status)

| # | Citation | Used for | Status |
|---|----------|----------|--------|
| 1 | Schild, Siegesmund, Vollbrecht & Mazurek (2001), *Characterization of granite matrix porosity and pore-space geometry by in situ and laboratory methods*, Geophys. J. Int. 146(1), 111–125, DOI 10.1046/j.0956-540x.2001.01427.x | intact granite porosity 0.5–1.2 % → lower bound / `phi_0` | VERIFIED (fetched) |
| 2 | Staněk & Géraud (2019), *Granite microporosity changes due to fracturing and alteration…*, Solid Earth 10, 251–274, DOI 10.5194/se-10-251-2019 | fractured/altered granite max measured porosity 6.5 % → "what 0.999 is not" | VERIFIED (fetched) |
| 3 | Heider (2021), *A review on phase-field modeling of hydraulic fracturing*, Eng. Fract. Mech. 253, 107881, DOI 10.1016/j.engfracmech.2021.107881 | phase-field framework; fully-cracked zone as fluid-filled void → upper-bound interpretation; matches the repo's permeability model citation | VERIFIED citation (fetched/searched); void-limit reading = modeling rationale |
| 4 | Westerly granite porosity ≈ 1.4 % (rock-deformation literature, multiple sources) | corroborating intact value | VERIFIED via search synthesis; corroborating only |
| 5 | Heider & Markert (2017), TPM phase-field hydraulic fracture in saturated porous media | supporting void/crack treatment | cited via ref. 3; UNVERIFIED at page level |

> Re-verify the DOIs in this table at finalization if the document is updated.
> Do not promote any `UNVERIFIED` row to an asserted fact without a resolving
> primary source.
