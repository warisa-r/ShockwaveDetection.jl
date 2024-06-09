# ShockDetection.jl

ShockDetection.jl is a Julia package for detecting shocks in time series data from Godunov solver.

## Features

- A method for shock detection.
- Easy-to-use functions for analyzing time series data.
- Compatible with Julia 1.9.

## Installation

You can install ShockDetection.jl using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode

# Functions

```@docs
read_output_file
convert_to_primitive
detect_shock
create_wave_animation
create_wave_animation_with_shock
```