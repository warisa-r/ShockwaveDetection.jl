using Documenter, ShockwaveDetection

makedocs(
    modules = [ShockwaveDetection],
    sitename = "ShockwaveDetection.jl Documentation",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ]
)