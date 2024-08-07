using Documenter
using UlamMethod

makedocs(
    sitename = "UlamMethod.jl",
    format = Documenter.HTML(),
    modules = [UlamMethod],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "API" => "api.md",
    ],
    warnonly = true
)

deploydocs(;
    repo = "github.com/70Gage70/UlamMethod.jl.git",
    target = "build", 
    versions = nothing
)
