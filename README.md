# ShockwaveDetection

[![Build Status](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/warisa-r/ShockwaveDetection.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Brief description

This package aims at processing the data from Euler's equation's solver [Euler2D.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) and utilizing [ShockwaveProperties.jl](https://github.com/STCE-at-RWTH/ShockwaveProperties.jl) to detect where the shock is.
Currently, it can:
- Detect 1D shock by detecting large gradients across properties with a customizable parameter `threshold`
- Detect shock curves in 2D and calculate normal vectors with various customizable parameters of the detection algorithms
- Visualize the change of properties along with shock positions
- Show performance metrics of each step of the detection

## Documentation

For detailed usage and API documentation, please refer to `doc/build/index.html`.

## Installation

To install this package, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.Registry.add("https://github.com/warisa-r/SWPRegistry.git")
Pkg.add(ShockwaveDetection)
