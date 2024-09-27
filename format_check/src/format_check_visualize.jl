using JuliaFormatter

function check_format(file::String)
    formatted_code = format_file(file)
    original_code = read(file, String)

    return formatted_code == original_code, formatted_code
end

function main()
    # Path to your visualize.jl file using a relative path
    file = joinpath(@__DIR__, "src", "visualize.jl")
    
    is_formatted, formatted_code = check_format(file)
    
    if !is_formatted
        println("Formatting issue found in: $file")
        println("Formatted code:")
        println(formatted_code)
        error("The file is not properly formatted.")
    else
        println("The file is properly formatted.")
    end
end

main()
