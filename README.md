This repository contains the mission-level analysis code, as well as the used debris dataset and the output dataset used in ["Feasibility Analysis of Small-Size Space Debris Removal in Low-Earth Orbit by Space-Based Laser Ablation"](https://www.researchgate.net/publication/380029470_Feasibility_Analysis_of_Small-Size_Space_Debris_Removal_in_Low-Earth_Orbit_by_Space-Based_Laser_Ablation) as presented at the 29th International Symposium on Space Flight Dynamics (ISSFD) 2024.

If you use any part of this repository or the accompanying paper, please cite us:
```
@conference{bogel_ldr,
  author       = {Elias Bögel and Hugo Buurmeijer and Lorenz Veithen and Frank Meijering and Gabriel Alves Teixeira and Daniel Rehling and Juan Bas Fernández and Per van Wolfswinkel and Niek Zandvliet and Jan Struzińksi},
  title        = {{Feasibility Analysis of Small-Size Space Debris Removal in Low-Earth Orbit by Space-Based Laser Ablation}},
  booktitle    = {International Symposium on Space Flight Dynamics (ISSFD)},
  year         = {2024},
  month        = {04},
  address      = {Darmstadt},
}
```


The two main studies done were by varying:
1) Range, ablation time and launch delay
2) Scan time, field of view and range

The configurations for both studies can be found in the the [configurations](https://github.com/eliasboegel/LDR/tree/main/configurations) directory.
Custom configurations files can be written in the same way as the provided configurations. Any parameter passed by the configuration overwrites the default value specified in [LDR.jl](https://github.com/eliasboegel/LDR/blob/main/LDR.jl).
The final time, as well as the final removal fraction is stored in [runs.csv](https://github.com/eliasboegel/LDR/blob/main/runs.csv), along with its input parameters.
This file contains the dataset used to produce figures 6 and 8 of the paper.
This data was generated on Julia 1.9.0.

# Instructions
1) Clone repository: `git@github.com:eliasboegel/LDR.git`
2) Run `packages.jl` to install the two required dependency packages
3) Run any of the configuration files
