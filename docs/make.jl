using Documenter
using EnergyStatistics

makedocs(
    sitename = "EnergyStatistics",
    format = Documenter.HTML(prettyurls = true),
    modules = [EnergyStatistics]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com:pfarndt/EnergyStatistics.jl.git"
)
