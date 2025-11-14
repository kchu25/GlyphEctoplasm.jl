"""
Distance key utilities for multi-motif analysis.
"""

"""
    distance_key_value(d)

Extract sortable tuple from distance key for lexicographic ordering of motif variants.
"""
function distance_key_value(d)
    if isa(d, Number)
        return (d,)
    elseif isa(d, Tuple)
        return d  # already a tuple
    elseif isa(d, NamedTuple)
        return Tuple(values(d))  # return all values as tuple for lexicographic sort
    elseif hasmethod(values, (typeof(d),))  # GroupKey and similar
        vals = collect(values(d))
        return isempty(vals) ? (0,) : Tuple(vals)  # return all values as tuple
    else
        try
            return (parse(Int, string(d)),)
        catch
            return (0,)
        end
    end
end

export distance_key_value
