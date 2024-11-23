using Documenter
using UlamMethod
import CairoMakie

makedocs(
    sitename = "UlamMethod.jl",
    format = Documenter.HTML(),
    modules = [UlamMethod],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "Making an App" => "app.md",
        "API" => "api.md"
    ],
    warnonly = true
)

deploydocs(;
    repo = "github.com/70Gage70/UlamMethod.jl.git",
    target = "build", 
    versions = nothing
)
