# PNG size optimization via median-cut color quantization.
#
# The plots this package emits (sequence logos, influence box/scatter, KDE
# panels) use only a handful of distinct colors, so re-encoding each truecolor
# PNG as an indexed (palette) PNG shrinks it by ~65% with no visible change.
# Crucially the filenames are unchanged, so no HTML/JS references need updating.
#
# Pure-Julia implementation (PNGFiles + IndirectArrays); no external binaries.

# Median-cut box over the *distinct* colors of the image (these plots have only
# a few thousand), each carrying a pixel count so splits/averages stay weighted
# by how much of the image each color covers.
mutable struct _ColorBox
    idxs::Vector{Int}     # indices into the unique-color list
    range::Int            # widest channel span (cached; -1 ⇒ needs recompute)
    ch::Int               # channel of that widest span
end

# Compute the count-weighted widest channel range of a box and cache it.
function _measure!(b::_ColorBox, colors::Vector{NTuple{3,UInt8}})
    if length(b.idxs) < 2
        b.range = -1; return
    end
    bestr, bestch = -1, 1
    @inbounds for ch in 1:3
        lo, hi = 0xff, 0x00
        for i in b.idxs
            v = colors[i][ch]; lo = min(lo, v); hi = max(hi, v)
        end
        r = Int(hi) - Int(lo)
        if r > bestr; bestr = r; bestch = ch; end
    end
    b.range = bestr; b.ch = bestch
end

# Median-cut over unique colors. `colors` are the distinct RGB tuples, `weights`
# their pixel counts. Returns (palette, palidx_of_color).
function _median_cut(colors::Vector{NTuple{3,UInt8}}, weights::Vector{Int}, ncolors::Int)
    root = _ColorBox(collect(1:length(colors)), 0, 1)
    _measure!(root, colors)
    boxes = _ColorBox[root]
    while length(boxes) < ncolors
        # pick the splittable box with the widest cached range
        best, bestrange = 0, -1
        for (bi, b) in enumerate(boxes)
            if b.range > bestrange; bestrange = b.range; best = bi; end
        end
        (best == 0 || bestrange < 0) && break    # nothing left to split
        b = boxes[best]
        ch = b.ch
        sort!(b.idxs, by = i -> colors[i][ch])
        # split at the weighted median so both halves carry ~equal pixel mass
        total = sum(@view weights[b.idxs])
        acc, cut = 0, 1
        @inbounds for (k, i) in enumerate(b.idxs)
            acc += weights[i]
            if acc * 2 >= total; cut = max(k, 1); break; end
        end
        cut = clamp(cut, 1, length(b.idxs) - 1)
        left = _ColorBox(b.idxs[1:cut], 0, 1)
        right = _ColorBox(b.idxs[cut+1:end], 0, 1)
        _measure!(left, colors); _measure!(right, colors)
        boxes[best] = left
        push!(boxes, right)
    end

    palette = Vector{NTuple{3,UInt8}}(undef, length(boxes))
    palidx_of_color = Vector{Int}(undef, length(colors))
    for (bi, b) in enumerate(boxes)
        r = g = bl = 0; w = 0
        @inbounds for i in b.idxs
            c = colors[i]; wi = weights[i]
            r += Int(c[1]) * wi; g += Int(c[2]) * wi; bl += Int(c[3]) * wi; w += wi
        end
        w = max(w, 1)
        palette[bi] = (round(UInt8, r / w), round(UInt8, g / w), round(UInt8, bl / w))
        @inbounds for i in b.idxs; palidx_of_color[i] = bi; end
    end
    return palette, palidx_of_color
end

# Flatten a loaded image to opaque RGB tuples, compositing any alpha over white
# (the page background) so transparency doesn't blow up the palette.
_rgb_tuple(c::Colorant) = begin
    rgba = RGBA{N0f8}(c)
    a = Float32(alpha(rgba))
    blend(x) = round(UInt8, 255 * (Float32(x) * a + (1 - a)))   # over white
    (blend(red(rgba)), blend(green(rgba)), blend(blue(rgba)))
end

"""
    quantize_png!(path; ncolors=64) -> Int

Re-encode the PNG at `path` in place as an indexed PNG with at most `ncolors`
palette entries (median-cut quantization). Returns the palette size used.

The image is overwritten only on success; on any failure the original file is
left untouched and a warning is emitted. Filenames are preserved, so existing
HTML/JS references keep working.
"""
function quantize_png!(path::AbstractString; ncolors::Int=64)
    try
        img = PNGFiles.load(path)
        h, w = size(img)
        flat = vec(img)

        # Collapse pixels to their distinct colors (+counts); quantize those.
        idmap = Dict{NTuple{3,UInt8},Int}()
        colors = NTuple{3,UInt8}[]
        weights = Int[]
        id_of_pixel = Vector{Int}(undef, length(flat))
        @inbounds for (k, c) in enumerate(flat)
            t = _rgb_tuple(c)
            id = get(idmap, t, 0)
            if id == 0
                push!(colors, t); push!(weights, 0)
                id = length(colors); idmap[t] = id
            end
            weights[id] += 1
            id_of_pixel[k] = id
        end

        palette, palidx_of_color = _median_cut(colors, weights, ncolors)
        index_of = [palidx_of_color[id] for id in id_of_pixel]
        pal = [RGB{N0f8}(reinterpret(N0f8, p[1]),
                         reinterpret(N0f8, p[2]),
                         reinterpret(N0f8, p[3])) for p in palette]
        ia = IndirectArray(reshape(index_of, h, w), pal)
        PNGFiles.save(path, ia)
        return length(pal)
    catch err
        @warn "quantize_png! skipped (kept original)" path exception=err
        return -1
    end
end

"""
    optimize_pngs!(dir; ncolors=64, verbose=true) -> NamedTuple

Recursively quantize every `*.png` under `dir` in place (see [`quantize_png!`]).
Returns `(; files, bytes_before, bytes_after)`. Safe to call more than once —
quantization is idempotent for already-indexed images.
"""
function optimize_pngs!(dir::AbstractString; ncolors::Int=64, verbose::Bool=true)
    files = 0
    before = 0
    after = 0
    for (root, _, names) in walkdir(dir)
        for name in names
            endswith(lowercase(name), ".png") || continue
            p = joinpath(root, name)
            b = filesize(p)
            quantize_png!(p; ncolors=ncolors) < 0 && continue
            files += 1
            before += b
            after += filesize(p)
        end
    end
    if verbose && before > 0
        pct = round(100 * (1 - after / before); digits=1)
        @info "optimize_pngs!: $files PNG(s) $(before >> 10) KB → $(after >> 10) KB (−$pct%)"
    end
    return (; files, bytes_before=before, bytes_after=after)
end
