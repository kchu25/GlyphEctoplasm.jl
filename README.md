# GlyphEctoplasm

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/GlyphEctoplasm.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/GlyphEctoplasm.jl/dev/)
[![Build Status](https://github.com/kchu25/GlyphEctoplasm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/GlyphEctoplasm.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/GlyphEctoplasm.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/GlyphEctoplasm.jl)

Motif visualization and interactive HTML rendering for sequence analysis.

## Overview

GlyphEctoplasm generates interactive HTML pages displaying sequence motifs with influence plots. It handles:

- Singleton motifs (single filter patterns)
- Multi-motifs (pairs, triplets, etc.)
- Mutation region analysis
- Banzhaf contribution visualization

## Installation

```julia
using Pkg
Pkg.add("GlyphEctoplasm")
```

## Usage

### Convolution-based motif analysis

```julia
using GlyphEctoplasm

config = ConvMotifConfig(data; filter_len=7, dpi=65, save_path="output")

json_motifs = init_json_dict()
html_dict = init_dict_for_html_render()

process_singletons!(df_singletons, config, json_motifs, html_dict)
process_multi_motifs!(dfs, config, json_motifs, html_dict; motif_size=2, group_id="pairs")

render_and_save_outputs!(json_motifs, html_dict, 1; save_path="output")
```

### Mutation region analysis

```julia
config = MutMotifConfig(data; filter_len=9, save_path="output")

plot_motifs_mut_case(data, model, df_filtered, dfs; save_path="output")
```

### High-level wrapper

```julia
plot_motifs_conv_case(data, model, motif_sizes, df_singletons, dfs, pts;
    save_path="results", page_title="My Analysis")
```

## Output

Generates an `index.html` with:
- Clickable motif cards showing sequence logos
- Influence distribution plots
- Modal dialogs with detailed views
- Navigation between motif groups
