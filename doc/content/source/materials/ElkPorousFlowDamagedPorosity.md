# ElkPorousFlowDamagedPorosity

!syntax description /Materials/ElkPorousFlowDamagedPorosity

## Description

`ElkPorousFlowDamagedPorosity` computes the damaged porosity stored in the property
`PorousFlow_porosity_qp_damaged` (or `PorousFlow_porosity_nodal_damaged` for a nodal
instance), used by the damage-coupled PorousFlow materials (e.g. the damaged Biot modulus).
Two update laws are selectable through `porosity_update_model`, and both are clamped to
`[porosity_lower_bound, porosity_upper_bound]`.

### `porosity_update_model = damage` (default)

The damage-driven, bounded-maximum porosity

\begin{equation}
\phi(d) = \phi_0 + (1 - \phi_0)\,\bigl[1 - (1 - d)^2\bigr],
\end{equation}

where $d$ is the phase-field damage ($d=0$ intact, $d=1$ fully broken). As $d \to 1$ the
porosity saturates at $1$ and is capped by `porosity_upper_bound` — the "bounded maximum
porosity" baseline. This reproduces the existing behavior; inputs that omit
`porosity_update_model` are unchanged.

### `porosity_update_model = strain`

The strain-based update of Liu et al. (2024, *Comput. Methods Appl. Mech. Engrg.* **429**,
117165), eq. (40):

\begin{equation}
\phi(\boldsymbol{\varepsilon}) = \phi_0 + \varepsilon_1,
\end{equation}

where $\varepsilon_1$ is the maximum (most-tensile) principal value of the kinematic strain
tensor `strain_property` (default `mechanical_strain`). It derives from a line crack of
aperture $\omega = h_e\,\varepsilon_1$ inside an element of edge $h_e$ (eqs. 32, 39), so the
cell-average porosity rises by $\omega/h_e = \varepsilon_1$. The update is **instantaneous
and reversible**: if the strain relaxes the porosity decreases again (unlike the monotone
damage law). In 2D plane strain $\varepsilon_1 \ge 0$, so pure compression returns
$\phi = \phi_0$.

#### Constraints for the strain law

- **Quadrature points only.** `mechanical_strain` is a qp property; selecting `strain` for a
  nodal porosity instance raises an error. Production inputs use the qp instance.
- **Requires a strain material.** A strain calculator (e.g. `ComputeSmallStrain`) must declare
  `strain_property`; otherwise the material errors at setup.
- **Non-AD, explicit coupling.** This material returns a plain `Real` with no strain
  derivative, so the strain-based porosity lags within a monolithic Newton iteration
  (consistent with how the damage-based porosity is already treated explicitly via the
  staggered phase-field sub-app).
- `mechanical_strain` is the kinematic (total) strain when no thermal/plastic eigenstrain is
  present. If such an eigenstrain is added, set `strain_property = total_strain` to keep
  eq. (40) representing the kinematic crack opening.

## Bounded-vs-strain comparison

To study the effect of the strain update against the bounded-maximum baseline, run the same
case twice changing only one parameter:

```
# baseline (bounded maximum porosity)
porosity_update_model = damage
# strain-based porosity (Liu 2024 eq. 40)
porosity_update_model = strain
```

and compare `PorousFlow_porosity_qp_damaged` and the porosity-dependent outputs (damaged
Biot modulus, fluid kinetic/elastic energy).

## Example Input Syntax

```
[Materials]
  [porosity_damaged]
    type = ElkPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
    porosity_update_model = strain   # eq. (40); omit or set 'damage' for the baseline
  []
[]
```

!syntax parameters /Materials/ElkPorousFlowDamagedPorosity

!syntax inputs /Materials/ElkPorousFlowDamagedPorosity

!syntax children /Materials/ElkPorousFlowDamagedPorosity
