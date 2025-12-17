
# hot accumulation helper
# accumulate_window! accumulates the window X[:, s:e, 1, didx] into `dest`.
# Accepts arbitrary Integer types for indices (e.g. Int32) and converts them to
# the platform `Int` before slicing. Returns the (possibly allocated) buffer
# which callers should reuse to avoid repeated allocations.
function accumulate_window!(dest::AbstractMatrix{T}, X, s::Integer, e::Integer, 
                           didx::Integer, buf::Union{Nothing,AbstractMatrix{T}}=nothing) where T<:AbstractFloat
    si = Int(s); ei = Int(e); didx_i = Int(didx)
    window = @view X[:, si:ei, 1, didx_i]

    if eltype(X) === T
        @inbounds dest .+= window
        return buf
    else
        if buf === nothing
            buf = Matrix{T}(undef, size(dest))  # undef is faster than similar
        end
        @inbounds @. buf = T(window)  # fused broadcast
        @inbounds @. dest += buf
        return buf
    end
end

# singleton wrapper (thin)
function build_singleton_count_matrices(gdf_filters::GroupedDataFrame, data, 
                                       filter_len::Int, ::Type{T}) where T<:AbstractFloat
    # Pre-allocate
    ngroups = length(gdf_filters)
    key_type = typeof(first(keys(gdf_filters)))
    res = Dict{key_type, Matrix{T}}()
    sizehint!(res, ngroups)
    
    @info "Building singleton count matrices for $(ngroups) filters..."
    @showprogress for k in keys(gdf_filters)
        g = gdf_filters[k]
        mat = zeros(T, 4, filter_len)
        buf = nothing
        
        # Indexed iteration
        @inbounds for i in axes(g, 1) # through rows
            s = g.position[i]
            e = s + filter_len - 1
            buf = accumulate_window!(mat, data.X, s, e, g.data_pt_index[i], buf)
        end
        res[k] = mat
    end
    return res
end

"""
build_count_matrices(gdf_d::GroupedDataFrame, data, motif_size, filter_len, ::Type{T}) -> Dict{Tuple{Vararg{Int}}, Matrix{T}}

Return a Dict mapping each distance-signature (immutable Tuple of Ints) to
the accumulated 4 x num_cols count matrix (eltype T). Uses a fast-path
when `eltype(data.X) === T` and otherwise converts into a preallocated buffer.
"""
function build_count_matrices_and_highlight(gdf_d::GroupedDataFrame, 
         data, motif_size::Int, 
         filter_len::Int, ::Type{T}) where T<:AbstractFloat

    # X = data.X
    X = data.onehot_sequences
    get_num_cols(mode) = sum(mode) + motif_size * filter_len
    
    # Pre-allocate dict with known size
    ngroups = length(gdf_d)
    key_type = typeof(first(keys(gdf_d)))
    key_to_mat = Dict{key_type, Matrix{T}}()
    sizehint!(key_to_mat, ngroups)

    # also make the highlighted regions
    highlighted_regions = Dict{key_type, Vector{UnitRange{Int}}}()
    for k in keys(gdf_d)
        highlighted_regions[k] = Vector{UnitRange{Int}}()
    end

    for k in keys(gdf_d)
        g = gdf_d[k]
        cols = get_num_cols(Tuple(Int(v) for v in values(k)))
        mat = zeros(T, 4, cols)
        buf = nothing
                
        # Use indexed iteration for clarity and potential performance
        @inbounds for i in 1:length(g.start)
            buf = accumulate_window!(mat, X, g.start[i], g.end[i], 
                                    g.data_pt_index[i], buf)
        end
        key_to_mat[k] = mat

        # highlighted regions
        highlighted_regions[k] = make_highlighted_region(motif_size, k, filter_len)
    end
 
    return key_to_mat, highlighted_regions
end


function make_highlighted_region(motif_size, d, filter_len::Integer)
    highlighted_region = Vector{UnitRange{Int}}()
    f_start = 1
    for idx = 1:motif_size
        f_end = f_start + filter_len - 1
        push!(highlighted_region, f_start:f_end)
        idx != motif_size && (f_start = f_end + d[idx] + 1)
    end
    return highlighted_region
end




### normalize count matrix
# mutating, type-stable, avoids temporaries; works with CPU/GPU arrays if you use similar(...)
function normalize_countmat!(out, countmat; pseudocount::Real=1e-3) 
    T = eltype(out)
    pc = T(pseudocount)

    @inbounds begin
        # convert + add pseudocount elementwise without allocating temporaries
        for i in eachindex(out)
            out[i] = T(countmat[i]) + pc
        end

        # column sums
        sums = sum(out, dims=1)
        # avoid divide-by-zero; use eps(T) or small value
        @. sums = max(sums, eps(T))

        # normalize in-place by broadcasting with the sums (fused on modern Julia)
        @. out = out / sums
    end

    return out
end

# convenience allocator wrapper
function normalize_countmat(countmat; pseudocount::Real=1e-3)
    # choose floating eltype if input is integer
    T = eltype(countmat) <: AbstractFloat ? eltype(countmat) : Float32
    out = similar(countmat, T)    # preserves device (CPU/GPU) and dims
    return normalize_countmat!(out, countmat; pseudocount=pseudocount)
end

# positional information saving
function save_pos_info_as_csv(df, filter_len, save_path)
    CSV.write(save_path, 
            DataFrame(seq_index  = df.data_pt_index,
                start_position = df.position,
                end_position   = df.position .+ filter_len .- 1,
                is_reversed_complement = zeros(Int32, length(df.position)) # for now
                )
    )
end

function save_pos_info_as_csv(flat_windows::Vector{NTuple{4,T}}, save_path) where T<:Integer
    n = length(flat_windows)
    df = DataFrame(seq_index = Vector{T}(undef, n),
                   start_position = Vector{T}(undef, n),
                   end_position = Vector{T}(undef, n),
                   is_reversed_complement = zeros(T, n)) 

    @inbounds for i in 1:n
        (didx, s, e, _) = flat_windows[i]
        df.seq_index[i] = didx
        df.start_position[i] = s
        df.end_position[i] = e
        df.is_reversed_complement[i] = 0 # for now
    end

    CSV.write(save_path, df)
end



"""
    save_as_meme(this_pfm, save_name::String; name = "", num_sites = 100)

Export position frequency matrix as MEME format file.

# Arguments
- `this_pfm`: Position frequency matrix (4 × width)
- `save_name`: Output file path
- `name`: Motif name (optional)
- `num_sites`: Number of sites parameter for MEME header

# MEME Format
```
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF name
letter-probability matrix: alength= 4 w= 8 nsites= 100 E= 0
 0.8 0.1 0.05 0.05
 0.1 0.7 0.1 0.1
 ...
```
"""
function save_as_meme(this_pfm, save_name::String; name = "", num_sites = 100)
    io = open(save_name, "w")
    print(io, "MEME version 4\n\n")
    print(io, "ALPHABET= ACGT\n\n")
    print(io, "strands: + -\n\n")
    print(io, "Background letter frequencies\n")
    print(io, "A 0.25 C 0.25 G 0.25 T 0.25\n\n")
    print(io, "MOTIF $name \n")
    print(
        io,
        "letter-probability matrix: alength= 4 w= $(size(this_pfm,2)) nsites= $num_sites E= 0\n",
    )
    for i in axes(this_pfm, 2)
        print(io, " ")
        for j = 1:4
            print(io, "$(this_pfm[j,i]) ")
        end
        print(io, "\n")
    end
    close(io)
end

# Build per-row vector-of-tuples: (data_pt_index, m_pos, m_pos+filter_len-1, 0)
# for each motif position column (m1_position, m2_position, ...).
# This fills a new column `:motif_windows` where each entry is a Vector of NTuple{4,Int}.
function build_motif_windows(d, motif_size::Integer, filter_len::Integer; 
    offset = 0)
    """
    Build motif windows for SubDataFrame `d`.

    Returns flat_windowns::Vector{NTuple{4,Int}} contains all windows consecutively,
    """
    pos_cols = m_position_symbols(motif_size)
    n = nrow(d)
    m = length(pos_cols)

    # snapshot columns for fast access
    pos_arrays = [d[!, c] for c in pos_cols]
    idx_col = d[!, :data_pt_index]

    total = n * m
    flat_windows = Vector{NTuple{4,Int}}(undef, total)
    p = 1
    @inbounds for i in 1:n
        didx = Int(idx_col[i])
        for j in 1:m
            pos = Int(pos_arrays[j][i]) + offset
            flat_windows[p] = (didx, pos, pos + filter_len - 1, 0)
            p += 1
        end
    end

    return flat_windows
end
