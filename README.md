# ShockwaveDetection

[![Build Status](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/warisa-r/ShockwaveDetection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/warisa-r/ShockwaveDetection.jl)


## Brief description

This package aims at processing the data from Euler's equation's solver [Euler2D.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and utilizing [ShockwaveProperties.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and detect where the shock is.
Currently, it can
- Detect 1D shock by detecting large gradients across properties with a customizable paramamer `threshold`
- Visualize the change of properties along with shock positions

This package also contains some documentation (far from complete or acceptable)

## Goals
- Develop 'better' methods to detect the shock
- Develop visualization of shock position in x axis against time
- Process 2D data (Waiting for the materials)
- Develop test series that alligns with the goal (probably from /examples)