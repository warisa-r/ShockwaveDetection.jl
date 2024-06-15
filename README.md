# ShockwaveDetection

[![Build Status](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Brief description

This package aims at processing the data from Euler's equation's solver [Euler2D.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and utilizing [ShockwaveProperties.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and detect where the shock is.
Currently, it can
- Detect 1D shock by detecting large gradients across properties with a customizable parameter `threshold`
- Visualize the change of properties along with shock positions

This package also an incomplete documentation (far from complete or acceptable)

## Goals
1. [x] Develop 'better' methods to detect the shock
   1. [ ] Make threshold a customizable fraction to be multiplied with a max gradient instead of just an unquantified number
2. [ ] Complete the Rankine-Hugoniot condition
    1. [x] Implement the rough check when uL and uR from state_behind are exactly the same position.
    2. [ ] Consider the smoothness by looking at CFL and grid size of hll solver 
3. [ ] Make a struct calling the solver to get the results or reading the output from a file and call ShockwaveDetectionAlgo with that
4. [ ] Develop visualization of shock position in the axis against time. Making velocity field a heat map? -> check GLMakie
5. [ ] Process 2D data (Waiting for the materials)
6. [ ] Develop a test series that aligns with the goal (probably from /examples)
