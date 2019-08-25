using Documenter, DocumenterMarkdown, JetPack

makedocs(
    format = Markdown(),
    sitename = "Foo"
)
cp("build/README.md", "../README.md", force=true)
