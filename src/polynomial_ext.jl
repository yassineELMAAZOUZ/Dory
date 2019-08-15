
export monomials_of_degree, dense_coefficients

@doc Markdown.doc"""
    monomials_of_degree(X::Vector{S} where S <: Hecke.MPolyElem{T} where T  ,itr) -->  
                                                                        ::Vector{S} where S <: Hecke.MPolyElem{T}

Given a list of variables X from a polynomial ring, return all monomials of degrees specified by `itr`.
`itr` can a list, iterator, or a single number.
"""
function monomials_of_degree(X::Vector{S} where S <: Hecke.MPolyElem{T} where T  ,itr)
    return vcat([monomials_of_degree(X,d) for d in itr]...)
end

function monomials_of_degree(X::Vector{S} where S <: Hecke.MPolyElem{T} where T ,d::Int64)

    # Edge case returns
    isempty(X) && error("Variables required in input")
    length(X) == 1 && return [X[1]^d]
    d < 0 && return Array{S,1}()
    d ==0 && return [parent(X[1])(1)]
    
    all_monomials = []
    n = length(X)
    Y = X[1:n-1]
    for j=0:d
        push!(all_monomials, [X[n]^j*m for m in monomials_of_degree(Y,d-j)])
    end
    return vcat(all_monomials...)
end

@doc Markdown.doc"""
    monomials_of_degree( R::Hecke.MPolyRing{T} where T,itr) --> ::Vector{S} where S <: Hecke.MPolyElem{T}

Given a polynomial ring, return all monomials of degrees specified by `itr`.
`itr` can a list, iterator, or a single number.
"""
function monomials_of_degree(R::Hecke.MPolyRing{T} where T, ds)
    X = gens(R)
    return monomials_of_degree(X,ds)
end


@doc Markdown.doc"""
    dense_coefficients(f::Hecke.MPolyElem{T}) where T <: RingElement

Return the list of monomial coefficients of the polynomial `f` for every monomial of degree up to
the total degree of `f`.
in the list L.
"""
function dense_coefficients(f::Hecke.MPolyElem{T}) where T <: RingElement
    R = parent(f)
    mons = monomials_of_degree(gens(R), 0:total_degree(f))
    return [coeff(f, m) for m in mons]
end


@doc Markdown.doc"""
    coeff(f :: Hecke.MPolyElem{T}, m :: Hecke.MPolyElem{T}) --> c :: T

Returns the coefficient of the monomial `m` in `f`. 

Note: docstring should be moved to Nemo.
"""
# function Hecke.coeff(f :: Hecke.MPolyElem{T}, m :: Hecke.MPolyElem{T}) where T
#     return Hecke.coeff(f, exponent_vector(m,1))
# end

@doc Markdown.doc"""
    coeff(f::Hecke.MPolyElem{T}, L) where T <: RingElement

Return the list of monomial coefficients of the polynomial `f` specified by the monomials
in the list L.
"""
function coeff(f::Hecke.MPolyElem{T}, L) where T <: RingElement
    return [Hecke.coeff(f,m) for m in L]
end
