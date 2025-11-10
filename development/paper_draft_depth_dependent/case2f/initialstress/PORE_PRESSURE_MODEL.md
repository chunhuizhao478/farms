# Pore Pressure Model - Summary

## Model Configuration

### Three Regions:

1. **Region 1 (0 to 6 km): Hydrostatic**
   - Pf = ρ_fluid · g · z = 1000 · 9.8 · z
   - Gradient: 9.8 kPa/m

2. **Region 2 (6 to 8 km): Transition Zone**
   - Linear interpolation from hydrostatic to overpressured
   - Smooth transition between Region 1 and Region 3

3. **Region 3 (> 8 km): Overpressured with constant gradient**
   - Pf = λ · ρ_rock · g · z = 0.9 · 2670 · 9.8 · z
   - Gradient: 23.55 kPa/m (90% of overburden gradient)
   - Same gradient as overburden, but scaled by λ

## Key Parameters:

- **lambda_pp = 0.9**: Pore pressure ratio (90% lithostatic)
- **A = 6000 m**: Transition start depth
- **B = 8000 m**: Transition end depth

## Physical Interpretation:

Below 8 km:
- **90% of overburden** supported by pore fluid pressure
- **10% of overburden** supported by rock skeleton (effective stress)
- Effective vertical stress: σ'_zz = -(1-λ) · ρ · g · z = -0.1 · ρ · g · z

## Example Values:

| Depth | Pf (MPa) | σ_zz (MPa) | σ'_zz (MPa) | % of overburden |
|-------|----------|------------|-------------|-----------------|
| 6 km  | 58.80    | -157.00    | -98.20      | 62.5%           |
| 7 km  | 123.60   | -183.16    | -59.56      | 32.5%           |
| 8 km  | 188.40   | -209.33    | -20.93      | 10.0%           |
| 10 km | 235.49   | -261.66    | -26.17      | 10.0%           |

## Adjusting the Model:

To change the effective stress level below 8 km, modify **lambda_pp**:

- lambda_pp = 1.0: Fully lithostatic (σ'_zz = 0)
- lambda_pp = 0.9: 10% effective stress (current)
- lambda_pp = 0.8: 20% effective stress
- lambda_pp = 0.5: 50% effective stress
- lambda_pp = 0.0: No pore pressure (σ'_zz = σ_zz)

## Output Plots:

The script generates two side-by-side plots:

1. **Left: Total Stress & Pore Pressure**
   - Shows |σ_zz|, |σ_xx|, |σ_yy|, σ_xy, and Pf
   - Pf gradient matches |σ_zz| gradient below 8 km

2. **Right: Effective Stress & Pore Pressure**
   - Shows effective stresses σ'_zz, σ'_xx, σ'_yy
   - Shows static and residual shear strengths
   - Pf shown as reference (dashed line)

