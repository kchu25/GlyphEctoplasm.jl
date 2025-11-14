using GlyphEctoplasm
using Documenter

DocMeta.setdocmeta!(GlyphEctoplasm, :DocTestSetup, :(using GlyphEctoplasm); recursive=true)

makedocs(;
    modules=[GlyphEctoplasm],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    sitename="GlyphEctoplasm.jl",
    format=Documenter.HTML(;
        canonical="https://kchu25.github.io/GlyphEctoplasm.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/GlyphEctoplasm.jl",
    devbranch="main",
)
