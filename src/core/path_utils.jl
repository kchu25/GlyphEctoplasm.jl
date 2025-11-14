# get_d_str(d) = join([i for i in d], "_")
# robust stringifier for group keys / mode keys
function get_k_mode_str(k)::String
    if k isa String
        return k
    end

    # helper that extracts the meaningful elements
    elems = try
        coll = collect(k)
        if !isempty(coll) && coll[1] isa Pair
            # DataFrames.GroupKey, NamedTuple iteration may yield Pair
            map(p -> p[2], coll)
        else
            coll
        end
    catch
        # fallback for scalars or non-iterables
        return string(k)
    end

    return join(map(string, elems), "_")
end

# Robust stringifier for distance / group keys
function get_d_str(x)::String
    try
        vals = collect(x)                      # works for Tuples, Vectors, NamedTuples (yields Pairs), GroupKey, etc.
        # If we got Pairs (key=>value), take the second element
        if !isempty(vals) && vals[1] isa Pair
            vals = map(v -> v[2], vals)
        end
        return join(map(string, vals), "_")
    catch
        # fallback for scalars or non-iterables
        return string(x)
    end
end

function get_descriptive_str(keys_, ds)
    ks = [(keys_[i], keys_[i+1]) for i = 1:length(keys_)-1]
    str = ""
    for ((m1, m2), d) in zip(ks, ds)
        str *= "pattern $(m1) and $(m2) are $(d) nucleotides apart<br>"
    end
    return str
end

"""
    get_filter_indices_str(k)

Extract filter indices from a NamedTuple key and format as a tuple string.
For singletons: returns just the number (e.g., "42")
For multi-motifs: returns tuple format (e.g., "(8, 112)" or "(8, 10, 112)")
"""
function get_filter_indices_str(k)
    try
        # Try to extract values from NamedTuple or similar
        vals = collect(values(k))
        if length(vals) == 1
            return string(vals[1])
        else
            return "(" * join(vals, ", ") * ")"
        end
    catch
        # Fallback for simple values
        return string(k)
    end
end





# get_pwm_save_path(k, d; t = 1) =
#     joinpath(get_pwm_save_folder(k; t = t), "$(get_d_str(d)).png")