# ShockwaveDetection

[![Build Status](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Brief description

This package aims at processing the data from Euler's equation's solver [Euler2D.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and utilizing [ShockwaveProperties.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and detect where the shock is.
Currently, it can
- Detect 1D shock by detecting large gradients across properties with a customizable parameter `threshold`
- Visualize the change of properties along with shock positions

This package also an incomplete documentation (far from complete or acceptable)

## Goals: Technical
1. [x] Develop 'better' methods to detect the shock
   1. [x] Make threshold a customizable fraction to be multiplied with a max gradient instead of just an unquantified number
2. [ ] Complete the Rankine-Hugoniot condition for 1D
3. [ ] Bug detected! In heatmap of velocity_field of 2D!!! Fix this!!!
4. [x] Develop visualization of shock position in the axis against time. Making velocity field a heat map? -> check GLMakie
5. [x] Process 2D data
6. [ ] Change the output format of shock positions to an x y coordinate the way [DBSCAN](https://github.com/JuliaStats/Clustering.jl) can take. Take a good look at [NearestNeighbor.jl](https://github.com/KristofferC/NearestNeighbors.jl)
7. What is wrong with the first frames, nothing is detected. Why is there a normal vector?
7. [ ] DBSCAN test
8. [ ] [LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl) seems to be a better choice due to more flexibility and more maintenance.
9. [ ] Develop a test series that aligns with the goal (probably from /examples)
   1. [ ] 1D: Check the consistency of the detection, does the shock positions disappear? What is the consistency rate? What is the velocity of the moving shock and is the detection of stationary shock consistent?

## Goals: Organizational
1. [ ] develop a formatting script
2. [ ] Update the documentation