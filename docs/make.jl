using Documenter
using UlamMethod

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "UlamMethod.jl",
    authors = "Gage Bonner",
    format = Documenter.HTML(),
    modules = [UlamMethod],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "Theory and Implementation" => "theory.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/70Gage70/UlamMethod.jl.git"
)
