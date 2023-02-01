using DocStringExtensions
using Documenter
using Polymers

makedocs(sitename = "Polymers", format = Documenter.HTML(), modules = [Polymers])

deploydocs(repo = "github.com/sandialabs/Polymers.git")
