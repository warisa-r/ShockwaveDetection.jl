using JuliaFormatter

# Provide the path to your clustering.jl file
file_to_check = "C:\\Users\\eivan\\Dropbox\\My PC (DESKTOP-AB2UFDR)\\Downloads\\ShockwaveDetection.jl-main (6)\\ShockwaveDetection.jl-main\\src\\cluster.jl"

@info "Checking format of $file_to_check"
try
    # Resolve the file path
    resolved_path = realpath(file_to_check)
    @info "Running format check on $resolved_path"

    # Format the file
    formatted_code = format(resolved_path)  # This will return the formatted code

    # Read the original code to check for differences
    original_code = read(resolved_path, String)

    # Check if the original and formatted codes are different
    if formatted_code == original_code
        @info "Formatting is good."
        println("Formatting is good.")
    else
        @warn "Formatting is bad. Please check the differences."
        println("Formatting is bad. Please check the differences.")
        println("\n--- Original Code ---\n")
        println(original_code)
        println("\n--- Formatted Code ---\n")
        println(formatted_code)
        
        # Optionally, you can overwrite the original file with formatted code
        # write(resolved_path, formatted_code)
    end
catch e
    @error "Error checking $file_to_check: $e"
    println("Error checking file: $e")
end
