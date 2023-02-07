using DocStringExtensions
using Documenter
using Polymers

makedocs(
    format = Documenter.HTML(),
    modules = [Polymers],
    pages = [
        "Home" => "index.md",
        "Modules" => ["Physics" => "physics.md"],
        "Examples" => [],
        "Indices and Tables" =>
            ["General Index" => "genindex.md", "Module Index" => "modindex.md"],
    ],
    sitename = "Polymers",
)

deploydocs(repo = "github.com/sandialabs/Polymers.git")
