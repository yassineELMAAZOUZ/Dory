push!(LOAD_PATH,"../src/")
using Documenter, DocumenterMarkdown, Dory, Pkg
makedocs(sitename="Dory.jl Documentation")

deploydocs(
  repo = "github.com/a-kulkarn/Dory.jl.git",
)

