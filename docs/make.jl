using Documenter
using UlamMethod

makedocs(
    sitename = "UlamMethod",
    format = Documenter.HTML(),
    modules = [UlamMethod]
)

deploydocs(;
    repo = "https://github.com/70Gage70/UlamMethod.jl"
)=
