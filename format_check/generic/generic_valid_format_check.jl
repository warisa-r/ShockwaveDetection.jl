using FilePaths

# Function to recursively check format of all .jl files in a directory
function check_format(directory::String)
    unformatted_files = String[]  # List to store unformatted files
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if endswith(file, ".jl")
                file_path = joinpath(root, file)
                # Check formatting
                result = run(`julia --check-bounds=no --color=yes --project -e "using JuliaFormatter; is_formatted(\"$file_path\")"`)
                if result != 0  # If the file is not formatted
                    push!(unformatted_files, file_path)
                end
            end
        end
    end

    return unformatted_files
end

# Specify the directory to search using relative path
# Assuming the script is located inside ShockwaveDetection.jl-main
project_directory = joinpath(@__DIR__, "src")  # Adjust "src" as needed

unformatted_files = check_format(project_directory)

if !isempty(unformatted_files)
    println("The following files are not formatted:")
    println(join(unformatted_files, "\n"))
    exit(1)  # Exit with an error code
end

println("All files are properly formatted.")
