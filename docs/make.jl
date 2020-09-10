using Documenter, JetPack

makedocs(sitename = "JetPack", modules=[JetPack])

deploydocs(
    repo = "github.com/ChevronETC/JetPack.jl.git",
)