"FARMS" : **F**ault **A**nd **R**upture **M**echanics **S**imulations

## Description

FARMS is a MOOSE-based application for simulating dynamic rupture and earthquake problem.

## Features

(Oct 1st 2023):

Various Frictional Laws: Slip Weakening, Rate-and-State Friction (including strong rate weakeing).

Complex Fault Geometry: Nonplanar Fault (main fault with branches / branch network).

Additional Constitutive Law: Continuum Damage-Breakage Rheology Model.

## Compile and Run

1. Install MOOSE: https://mooseframework.inl.gov/getting_started/installation/ 

2. Inside ``` ~/projects/ ``` folder, clone the repository: ``` git clone https://github.com/chunhuizhao478/farms.git ```

3. Inside ``` ~/projects/farms ``` folder, run the code: ``` make -j8 ```

4. To run code, for example: 

For single core: ``` ./farms-opt -i examples/benchmark_tpv2052D/tpv2052D.i ```

For multiple cores: ``` mpirun -np #[num of cores] ./farms-opt -i examples/benchmark_tpv2052D/tpv2052D.i ```

## Contributors

Chunhui Zhao and Mohamed Abdelmeguid, Ahmed Elbanna, PhD

Mechanics of Complex Systems Group\
Department of Civil and Environmental Engineering\
University of Illinois Urbana-Champaign

## Contact Information

Chunhui Zhao (chunhui3@illinois.edu)
