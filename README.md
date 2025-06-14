"FARMS" : **F**ault **A**nd **R**upture **M**echanics **S**imulations

## Description

FARMS is a MOOSE-based application (https://github.com/idaholab/moose.git) for simulating dynamic rupture and earthquake problem. 

## Features

- Various Frictional Laws: Slip Weakening, Rate-and-State Friction with SCEC Benchmark Verification.

![TPV2053DVEL](https://github.com/user-attachments/assets/78e71f62-7a27-4783-bf20-9f87214500d8)

- Complex Fault Geometry: Nonplanar Fault (main fault with branches / branch network).

![2023 Turkey-Syria Earthquake](https://github.com/user-attachments/assets/23b885e1-e0a7-4c1d-9ac6-dff92372ba25)

- Kinematics: Small-strain and Finite deformation

- Additional Constitutive Law: Continuum Damage-Breakage Rheology Model.

![Damage-Breakage Supershear](https://github.com/user-attachments/assets/0cb608fb-0dfd-43ab-bc76-5bc4eb4ce3e6)

- Multphysics Coupling: Poroelasticity

## Compile and Run

1. Install MOOSE: https://mooseframework.inl.gov/getting_started/installation/ 

2. Inside ``` ~/projects/ ``` folder, clone the repository: ``` git clone https://github.com/chunhuizhao478/farms.git ```

3. Inside ``` ~/projects/farms ``` folder, run the code: ``` make -j8 ```

4. To run code, for example: 

For single core: ``` ./farms-opt -i examples/benchmark_tpv2052D/tpv2052D.i ```

For multiple cores: ``` mpirun -np #[num of cores] ./farms-opt -i examples/benchmark_tpv2052D/tpv2052D.i ```

## Contributors

Chunhui Zhao, Mohamed Abdelmeguid, Amr Ibrahim and Ahmed Elbanna, PhD

Mechanics of Complex Systems Group\
Department of Civil and Environmental Engineering\
University of Illinois Urbana-Champaign

## Contact Information

Chunhui Zhao (chunhui3@illinois.edu)
