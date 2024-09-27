using JuliaFormatter

# Provide the relative path to your fitting.jl file
relative_path_to_check = joinpath("src", "fitting.jl")

# Get the current working directory
current_dir = pwd()

# Combine current directory and relative path
file_to_check = joinpath(current_dir, relative_path_to_check)

@info "Checking format of $file_to_check"
try
    # Resolve the file path
    resolved_path = realpath(file_to_check)
    @info "Running format check on $resolved_path"
    format(resolved_path)  # Run the format check on the resolved path
catch e
    @error "Error checking $file_to_check: $e"
end
