using Documenter
using EnergyStatistics

makedocs(
    sitename = "EnergyStatistics",
    format = Documenter.HTML(),
    modules = [EnergyStatistics]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
