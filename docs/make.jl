push!(LOAD_PATH,"../src/")
using Pkg, Documenter, DocumenterMarkdown, Dory
makedocs(sitename="Dory.jl Documentation")

deploydocs(
  repo = "github.com/a-kulkarn/Dory.git",
)

