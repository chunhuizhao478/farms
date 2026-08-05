# MOOSE syntax compatibility check

The published deck was written against MOOSE `9de1ffd476` (2023-11-09). This branch pins MOOSE as
a submodule at `8fd9652711e74258b48d799bc8819bfa9bb4b33b`, merged **2 July 2024**
(PR #28035, "Add an option to not contribute to RC coefficients").

Three constructs in the deck were deprecated between those two versions. Each was checked against
the pinned commit before deciding whether the deck needed rewriting.

| Construct in the deck | Status at the pinned MOOSE | Action taken |
| --- | --- | --- |
| `Modules/TensorMechanics/Master` | Registered via `registerDeprecatedSyntax("CommonSolidMechanicsAction", "Modules/TensorMechanics/Master", ...)` in `modules/solid_mechanics/src/base/SolidMechanicsApp.C`. Works, warns. | **Kept as published.** |
| `Modules/TensorMechanics/CohesiveZoneMaster` | Registered via `registerDeprecatedSyntax("CommonCohesiveZoneAction", "Modules/TensorMechanics/CohesiveZoneMaster", ...)` in the same file. Works, warns. | **Kept as published.** |
| `Outputs/interval` | `params.addParam<unsigned int>("interval", ...)` plus `params.deprecateParam("interval", "time_step_interval", "02/01/2025")` in `framework/src/outputs/Output.C`. Works, warns. | **Kept as published.** |

**Conclusion: no syntax migration is required.** The deck runs on the pinned MOOSE unmodified.
This was the main open risk when this example was planned, and it resolved in favour of shipping
the published deck verbatim.

## Expected warnings — all benign

A correct run prints these. They are not errors and do not affect results:

```
*** Warning, This code is deprecated and will be removed in future versions:
  .../test_planarfault_main.i:256.5: The 'Modules/TensorMechanics/CohesiveZoneMaster' syntax is
  deprecated. Please use 'Physics/SolidMechanics/CohesiveZone' instead.

*** Warning, This code is deprecated and will be removed in future versions:
  .../test_planarfault_main.i:266.7: The 'Modules/TensorMechanics/Master' syntax is deprecated.
  Please use 'Physics/SolidMechanics/QuasiStatic' instead.

*** Warning, This code is deprecated and will be removed in future versions:
  .../test_planarfault_main.i:574: (Outputs/interval): 'interval' has been deprecated and will be
  removed on 02/01/2025. Please use 'time_step_interval' instead.
```

You will also see, from the solid-mechanics action:

```
*** Info *** SolidMechanics Action: selecting 'total small strain' formulation.
```

## If you bump the MOOSE submodule

The `interval` parameter carries a stated removal date of **2025-02-01**, which the pinned
submodule predates. If you advance the submodule past that point, this deck will stop parsing.

The fix is a pure rename with no effect on results:

```
# before
[Outputs]
  exodus = true
  interval = 200

# after
[Outputs]
  exodus = true
  time_step_interval = 200
```

The two `Modules/TensorMechanics/*` blocks will need the corresponding
`Physics/SolidMechanics/{QuasiStatic,CohesiveZone}` migration, which is **not** a pure rename —
verify against `reference/fingerprint_alpha_B.csv` after any such change, because the newer
action defaults differ.

Do not make either change speculatively. The deck as committed matches the published run, and
that is the point of this folder.
