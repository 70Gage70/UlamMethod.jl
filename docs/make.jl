using Documenter
using UlamMethod

makedocs(
    sitename = "UlamMethod.jl",
    authors = "Gage Bonner",
    format = Documenter.HTML(),
    modules = [UlamMethod],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/70Gage70/UlamMethod.jl.git",
    versions = nothing
)
