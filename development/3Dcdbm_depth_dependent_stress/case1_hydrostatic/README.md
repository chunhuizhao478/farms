Case 1: Hydrostatic Pore Pressure
Created By Chunhui Zhao, Nov 11th, 2025

- To check initial stress components distribution: python3 development/3Dcdbm_depth_dependent_stress/case1_hydrostatic/initialstress/plotsts.py

- Static solve:
(1) "InitialStressStrainTPV26" implements the stress components as shown in plotsts.py
(2) The problem setup includes several parts:
1. boundary tractions: see [BCs] section
2. eigenstrain: convert stress into strain using "ComputeEigenstrainFromInitialStress" and add it to "ComputeSmallStrain"
3. for overbuden stress, we use body force + fixed bottom box z direction to maintain equilibrium
4. manually check in paraview if the imposed stress (static_initial_stress_tensor_xx) is equal (or close to) resolved stress (stress_xx)

- Dynamic solve:
(1) The nucleation is done by reducing the fault strength, see SCEC TPV26
(2) local xi -- ElkRadialAverageUpdated(nonlocal xi) -- ElkNonlocalEqstrainUpdated(extract) -- ComputeDamageBreakageStress3DSlipWeakeningNonlocal (called)
