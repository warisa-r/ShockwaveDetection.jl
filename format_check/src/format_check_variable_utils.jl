using JuliaFormatter

# Get the path of the current script (or use the current working directory)
current_dir = dirname(@__FILE__)  # If running as a script file, use @__FILE__
# For running directly in REPL, use `pwd()` to get the current working directory
# current_dir = pwd()

# Construct the relative path to the file based on the current directory
file_path = joinpath(current_dir, "src", "variable_utils.jl")

# Read the original file content
original_content = read(file_path, String)

# Format the content without changing the file
formatted_content = format_text(original_content, margin=92)  # Adjust margin if needed

# Check if the original content matches the formatted content
if original_content == formatted_content
    println("The file is properly formatted.")
else
    println("The file is not properly formatted.")
end
