using Documenter, DocumenterVitepress
using UlamMethod

makedocs(
    sitename = "UlamMethod.jl",
    authors = "Gage Bonner",
    format = MarkdownVitepress(
        repo = "github.com/70Gage70/UlamMethod.jl.git",
    ),
    modules = [UlamMethod],
    pages = [
        "Home" => "index.md",
        "Advanced Usage" => "advanced.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/70Gage70/UlamMethod.jl.git",
    target = "build", # this is where Vitepress stores its output
    versions = nothing
)
