
@doc Markdown.doc"""
    maximum(abs, A::SMat{fmpz}) -> fmpz
  Finds the largest, in absolute value, entry of $A$. 
  Note that `abs` must be defined for the entry type.
"""

import Base.maximum
function maximum(::typeof(abs), A::SMat{T} where T)
    if length(A.rows) == 0
        return zero(FlintZZ)
    end
    m = abs(A.rows[1].values[1])
    for i in A.rows
        for j in i.values
            if cmpabs(m, j) < 0
                m = j
            end
        end
    end
    return abs(m)
end


function size(A::SMat{T} where T, i::Int64)
    return size(A)[i]
end

function rows(A::SMat{T} where T)
    return A.rows
end

## Temporary patch.
function Base.deepcopy(A::SMat{T} where T)
    return Hecke.copy(A)
end
