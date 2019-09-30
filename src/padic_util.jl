

# simple function to invert a permutation array.
function inverse_permutation(A::Array{Int64,1})
    Pinv = fill(0,length(A))
    for i=1:length(A)
        Pinv[A[i]] = i
    end
    return Pinv
end

##############################################################################################
#                                                                                            #
#                          Basic extension to padic numbers                                  #
#                                                                                            #
##############################################################################################

import Base: +, abs 
import Hecke.valuation

function +(x::padic) return x end
function /(x::padic,y::padic) return x//y end

# Access to the precision fields.
"""

"""
function precision(Qp::FlintPadicField) return Qp.prec_max end
function prec(Qp::FlintPadicField) return Qp.prec_max end

    
# Potential fix for the bug with the valuation function
# Note: there may be issues if Hecke depends on the
# valuation being non-infinite.
#
# function valuation(x::padic)
#     if iszero(x)
#         return Inf
#     end
#     return Int64(x.v)
# end

# typesafe version
function float64_valuation(x::padic)
    if iszero(x)
        return Inf
    end
    return Float64(x.v)
end


function abs(x::padic)
    p = parent(x).p
    return Float64(p)^(-valuation(x))
end

function modp(x::padic)
    Fp = ResidueRing(FlintZZ,parent(x).p)
    return Fp(lift(x))
end

## Test utilities
function test_rings()
    return Qp = FlintPadicField(7,20), ResidueRing(FlintZZ,7)
end

import Base.rand
function rand(Qp::FlintPadicField)
    p = Qp.p
    N = Qp.prec_max
    return Qp(rand(1:BigInt(p)^N))*Qp(p)^(rand(-N:N))
end

function rand_padic_int(Qp::FlintPadicField)
    p = Qp.p
    N = Qp.prec_max
    return Qp(rand(1:BigInt(p)^N))
end


function random_test_matrix(Qp,n=4)
    A = matrix(Qp, fill(zero(Qp),n,n))
    for i=1:n
        for j=1:n
            A[i,j] = rand_padic_int(Qp)
        end
    end
    return A
end

##############################################################################################
#                                                                                            #
#                          Polynomials over p-adic fields                                    #
#                                                                                            #
##############################################################################################

# Getting coefficients of a flint polynomial is not intuitive.
# Flint system crashes if coefficients are not integers.
# Flint system crashes if leading coefficients are divisible by p.

# Lift termwise to a polynomial over the Flintegers.
import Hecke.lift
function lift(f :: Hecke.Generic.Poly{padic})
    R,_ = PolynomialRing(FlintZZ)
    return R([lift(c) for c in f.coeffs])
end

# This function is...kind of a hack.
# It is also very buggy since FLINT can only handle a specific case
# (integer polynomial, non-vanishing leading coefficient mod p)
function factor(f :: Hecke.Generic.Poly{padic})
    QpX = f.parent
    Qp = QpX.base_ring
    N = prec(Qp)
    
    f_int = lift(f)
    H = factor_mod_pk_init(f_int,Qp.p)
    D = factor_mod_pk(H,N)

    return Dict( QpX(lift(k))=>D[k] for k in keys(D))   
end

##############################################################################################
#                                                                                            #
#                               p-adic linear algebra                                        #
#                                                                                            #
##############################################################################################

# Compute a factorization of a padic matrix A = QR, where R is
# upper triangular (Borel) and Q lies in the maximally compact
# subgroup of SL(Qp) = SL(Zp).
#
# It turns out that the algorithm to do this is just the LU factorization
# with pivoting.

struct QRPadicPivoted
    Q::Hecke.Generic.MatElem{padic}
    R::Hecke.Generic.MatElem{padic}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

struct QRPadicSparsePivoted
    Q::Hecke.SMat{padic}
    R::Hecke.SMat{padic}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

@doc Markdown.doc"""
    padic_qr(A :: Hecke.Generic.MatElem{padic} ; col_pivot :: Union{Val{true},Val{false}}) --> F :: QRPadicPivoted

The return type of `F` is a QRPadicPivoted, with fields `F.Q, F.R, F.p, F.q` described below.
                 
Compute the p-adic QR factorization of `A`. More precisely, compute matrices `Q`,`R`, and an arrays `p`, `q` such that 

    A[F.p,F.q] = Q*R

If `col_pivot=Val(false)`, then `F.q = [1,2,...,size(A,2)]`.

#-------------------

INPUTS:
`A`         -- a matrix over Qp
`col_pivot` -- a type, either `Val(true)` or `Val(false)`, indicating whether column permutations
             should be used to move p-adically large entries to the pivot position.

"""
function padic_qr(A::Hecke.Generic.MatElem{padic};
                  col_pivot=Val(false) :: Union{Val{true},Val{false}})

    # Set constants
    n = size(A,1)::Int64
    m = size(A,2)::Int64
    basezero = zero(A.base_ring)
    
    L= identity_matrix(A.base_ring,n)
    Lent = L.entries::Array{padic,2}  
    Umat= deepcopy(A)
    U = Umat.entries
    
    P= Array(1:n)
    Pcol=Array(1:m)

    # We cache the maximum value of the matrix at each step, so we save an iteration pass
    # through the matrix.
    val_list = float64_valuation.(U)
    min_val, min_val_index = findmin( val_list );
    
    # Allocate specific working memory for multiplications.
    container_for_swap    = padic(U[1,1].N)
    container_for_product = padic(U[1,1].N)
    container_for_div     = padic(U[1,1].N)
    
    # Allocate a 2-element array to hold the index of the maximum valuation.
    min_val_index_mut = [x for x in min_val_index.I]

    # Allocate a function to zero consecutive entries of a column
    function zero_subdiagonal_of_column!(U,k::Int64)
        for j = k+1:n
            zero!(U[j,k])
        end
    end
    
    
    for k=1:(min(n,m)::Int64)

        if col_pivot==Val(true)
            col_index=min_val_index_mut[2]
            if col_index!=k
                # interchange columns m and k in U
                for r=1:n
                    U[r,k], U[r,col_index] = U[r,col_index], U[r,k]
                end
                
                # interchange entries m and k in Pcol
                Pcol[k], Pcol[col_index] = Pcol[col_index], Pcol[k]
            end
        end
        
        val_list = float64_valuation.(U[k:n,k])
        minn, row_pivot_index = findmin( val_list );
        if minn==Inf continue end

        row_pivot_index=row_pivot_index+k-1;
        if row_pivot_index!=k

            # interchange rows `row_pivot_index` and `k` in U
            for r=1:m
                U[k,r], U[row_pivot_index,r] = U[row_pivot_index,r], U[k,r]
            end               
            
            # interchange entries `row_pivot_index` and k in P
            P[k],P[row_pivot_index] = P[row_pivot_index],P[k]

            # swap columns corresponding to the row operations already done.
            swap_prefix_of_row!(Lent, k, row_pivot_index)
        end

        # Reset min_valuation for selecting new pivot.
        min_val = Inf

        # Note to self: with Julia, the optimal thing to do is split up the row operations
        # and write a for loop.
        # The entries left of the k-th column are zero, so skip these.
        # Cache the values of L[j,k] first.
        #
        if iszero(U[k,k])
            # If col_pivot == true, then we don't need to perform further column swaps
            # in this case, since the element of largest valuation in the lower-right
            # block is zero. In fact, no further operations need to be performed.
            continue
        end 

        # The use of the inversion command preserves relative precision. By row-pivoting,
        # the extra powers of p cancel to give the correct leading term.
        # the "lost" digits of precision for L[j,k] can simply be set to 0.
        container_for_inv = inv(U[k,k]) 
        
        for j=k+1:n
            Hecke.mul!(L[j,k],U[j,k], container_for_inv)
            L[j,k].N = parent(L[j,k]).prec_max            # L[j,k] is really an integer.
        end

        zero_subdiagonal_of_column!(U,k)
        
        for r=k+1:m
            for j=k+1:n
                # Compute U[j,r] = U[j,r] - L[j,k]*U[k,r]                
                Hecke.mul!(container_for_product, L[j,k], U[k,r])
                _unsafe_minus!(U[j,r], container_for_product)
                
                # Update the smallest valuation element
                if float64_valuation(U[j,r]) < min_val
                    min_val = float64_valuation(U[j,r])
                    min_val_index_mut[1] = j
                    min_val_index_mut[2] = r
                end
            end
        end
    end

    @assert iszero(A[P,Pcol] - L*Umat)
    
    return QRPadicPivoted(L,Umat,P,Pcol)
end

function swap_prefix_of_row!(Lent, k::Int64, i::Int64)
    # The index of the diagonal point is (k,k)
    for r=1:(k-1)
        container_for_swap = Lent[k,r]
        Lent[k,r] = Lent[i,r] 
        Lent[i,r] = container_for_swap
    end
    return
end

# Performs subtraction in-place, x-> x-y 
function _unsafe_minus!(x::padic, y::padic)
    x.N = min(x.N, y.N)
    ccall((:padic_sub, :libflint), Nothing,
          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
          x, x, y, parent(x))
    return
end

#######################################################################################

# Sparse QR-algorithm

#######################################################################################

### "Seriously!? how is this not implemented?" block.
import Base.findfirst
function Base.findfirst(A::Array{Int64,1}, j::Int64)
    return findfirst(x->(x==j), A)
end

@doc Markdown.doc"""
    setindex!(A::SMat{T}, a::T, i::Int64, j::Int64) -> nothing
Given a sparse matrix $A$, set $A[i,j] = a$.
"""
function Base.setindex!(A::SMat{T}, a::T, i::Int64, j::Int64) where {T<:Hecke.RingElem}
    i < 1 && error("Index must be positive")
    srow = A[i]
    
    p = findfirst(isequal(j), srow.pos)
    if p === nothing
        irange = searchsorted(srow.pos,j)
        splice!(srow.pos, irange, [j])
        splice!(srow.values, irange, [a])
        A.nnz += 1
    else
        srow.values[p] = a
    end
    A[i] = srow
    return
end

@doc Markdown.doc"""
    sparse_identity(R::T, n::Int64) where T<:Hecke.Ring --> Id

Given a ring `R`, and an integer `n`, return the `n x n`-identity matrix over the ring R.

#-------------------

INPUTS:
R         -- type ::Hecke.Ring
col_pivot -- a type, either Val(true) or Val(false), indicating whether column permutations
             should be used to move p-adically large entries to the pivot position.

"""
function sparse_identity(R::T, n::Int64) where T<:Hecke.Ring
    II = Hecke.sparse_matrix(R)
    for i=1:n
        srow = sparse_row( R, Array{Int64,1}([i]), Array{padic,1}([R(1)]) )
        push!(II, srow)
    end
    return II
end
### End block of misery.

# The index of the diagonal point is (k,k)
function swap_prefix_of_column!(L, diagonal_index::Int64, i::Int64)
    k = diagonal_index
    for r=1:k-1
        L[r,k],L[r,i] = L[r,i], L[r,k]
    end
    return
end


"""
    padic_qr(A; col_pivot) --> F

The return type of `F` is a QRPadicPivoted, with fields `F.Q, F.R, F.p, F.q` described below.
                 
Compute the p-adic QR factorization of A. More precisely, compute matrices `Q`,`R`, and an arrays `p`, `q` such that 

    A[F.p,F.q] = Q*R

If col_pivot=Val(false), then F.q = [1,2,...,size(A,2)].

#-------------------

INPUTS:
A         -- a matrix over Qp, A::Hecke.Generic.MatElem{padic}
col_pivot -- a type, either Val(true) or Val(false), indicating whether column permutations
             should be used to move p-adically large entries to the pivot position.

"""
function padic_qr(A::Hecke.SMat{padic};
                  col_pivot=Val(false) :: Union{Val{true},Val{false}})
    
    # Set constants
    Qp = A.base_ring
    n = size(A,1)::Int64
    m = size(A,2)::Int64
    
    # We store the ***transpose*** of L as a sparse matrix, and flip everything at the end.
    # Allocate the rows of Ltrans ahead of time.
    Ltrans = sparse_identity(Qp, n)
    
    U= deepcopy(A)    
    P= Array(1:n)
    Pcol=Array(1:m)

    # Function to pivot and return the list of rows with a non-zero entry at index k.
    # in the subdiagonal.
    function pivot_and_select_row_indices(U, Ltrans, k, piv)

        # Scan through the matrix to check if a column swap is needed.
        if col_pivot==Val(true)
            minn = Inf
            mindex = piv
            for j=k:n
                srow = U[j]
                if isempty(srow) break end
                
                rowmin, rowmindex = findmin( valuation.(srow.values) )
                if rowmin < minn
                    minn = rowmin
                    mindex = srow.pos[rowmindex]
                end
            end

            if mindex != piv
                Hecke.swap_cols!(U,piv,mindex)
                Pcol[piv], Pcol[mindex] = Pcol[mindex], Pcol[piv]
            end
        end

        # Scan through the matrix to check if a rowswap is needed.
        valuation_index_pairs = [ (valuation(U[j, piv]), j) for j=k:n if !iszero(U[j, piv]) ]        
            
        if !isempty(valuation_index_pairs)

            row_pivot_index = last(minimum(valuation_index_pairs))
            rows_with_entry_at_piv = last.(valuation_index_pairs)
            
            if row_pivot_index!=k
                Hecke.swap_rows!(U,k,row_pivot_index)
                P[k],P[row_pivot_index] = P[row_pivot_index],P[k]

                # swap columns corresponding to the row operations already done.
                # Do not swap the diagonal elements.
                swap_prefix_of_column!(Ltrans, k, row_pivot_index)
            end

            # Remove k or leading zeros from the list of rows to iterate over.
            # (Leading zeros can be introduced by a swap.)
            filter!(j->(j!=k && !iszero(U[j,piv])), rows_with_entry_at_piv)
        else        
            return last.(valuation_index_pairs)
        end
    end ### END PIVOT FUNCTION

    # Initialize the shift index to determine the critial pivot location.
    shift=0
    k=1
    
    @time while k <= (min(n,m)::Int64) && k+shift <= m::Int64

        ### Pivot logic ###
        
        # set the pivot column and determine the rows to apply elimination.
        piv = k + shift
        rows_with_entry_at_piv = pivot_and_select_row_indices(U, Ltrans, k, piv )

        # If everything is zero, shift the algorithm to operate on the
        # right rectangular submatrix window.
        if isempty(rows_with_entry_at_piv) shift+=1; continue; end   

        
        ### Row operations loop ###
        container_for_inv = inv(U[k,piv])
        
        for j in rows_with_entry_at_piv
                        
            # The "lost" digits of precision for L[j,k] can simply be set to 0.
            # as L[j,k] is really an integer.
            Ltrans[k,j] = U[j,piv]*container_for_inv

            if Ltrans[k,j].N < prec(Qp)
                Ltrans[k,j].N = prec(Qp)
            end
            
            if Ltrans[k,j] != 0
                Hecke.add_scaled_row!(U, k, j, -Ltrans[k,j])
            elseif valuation(Ltrans[k,j]) < prec(Qp)
                error("Problem related to Hecke's `add_scaled_row` function encountered.")
            end
        end

        #Update loop counter.
        k += 1        
    end

    

    # The test is actually quite expensive, but we keep it for now.
    @time @assert iszero( matrix(A)[P,Pcol] - transpose(matrix(Ltrans))*matrix(U) )

    return QRPadicSparsePivoted( transpose(Ltrans),U,P,Pcol)
end


########################################################################################

# IMPORTANT!
# We deviate slightly from LinearAlgebra's SVD structure by putting a diagonal matrix for S.
struct SVDPadic
    U::Hecke.Generic.MatElem{padic}
    S::Hecke.Generic.MatElem{padic}
    Vt::Hecke.Generic.MatElem{padic}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

# A padic analogue for svd
@doc Markdown.doc"""
    svd(A :: Hecke.Generic.Matelem{padic}) -> SVDPadic

  Compute the singular value decomposition (SVD) of A and return an SVDPadic object.

  A pAdic singular value decomposition is a factorization of the form

  A = U * S * Vt

  where U,Vt are matrices in GL_n(Zp). For efficiency reasons, we give a factorization
  as well as two arrays `F.p` and `F.q` such that

  A[p, q] = U*S*Vt

  where `U` is lower triangular with ones on the diagonal and `Vt` is upper triangular
  with ones on the diagonal. The singular values in S are sorted in descending order of padic 
  absoute value.

  The return types are
     U, S, Vt :: Hecke.Generic.MatElem{padic}
     p,q      :: Array{Int64,1}

"""
function svd(A::Hecke.Generic.MatElem{padic})

    F = padic_qr(A, col_pivot=Val(true))
    G = padic_qr(transpose(F.R))

    @assert G.p == [i for i=1:length(G.p)]

    U = deepcopy(F.Q)
    S = transpose(G.R)
    Vt= transpose(G.Q)
    
    @assert iszero( A[F.p,F.q] - U*S*Vt)

    return SVDPadic(U,S,Vt,F.p,F.q)
end

# stable version of rank for padic matrices.
@doc Markdown.doc"""
    rank(A::Hecke.Generic.MatElem{padic})

  Compute the rank of a padic matrix by counting how many singular values satisfy `iszero(a)`.
"""
function rank(A::Hecke.MatElem{padic})
    n = nrows(A)
    m = ncols(A)
    F = padic_qr(A)

    rank=0
    for i=1:min(n,m)
        if !iszero(F.R[i,i])
            rank += 1
        end
    end
    return rank
end

# Returns the p-adic singular values of a matrix
@doc Markdown.doc"""
    singular_values(A::Hecke.MatElem{padic}) -> Array{padic, 1}

Returns the list of diagonal elements in the singular value decomposition of the matrix `A`.
"""
function singular_values(A::Hecke.MatElem{padic})
    F = padic_qr(A,col_pivot=Val(true))
    return [ F.R[i,i] for i=1:minimum(size(A)) ]
end

# stable version of nullspace for padic matrices.
import Hecke.nullspace
@doc Markdown.doc"""
    nullspace(A::Hecke.MatElem{padic}) -> (nu, N)

Computes the nullspace of a padic matrix `A`. The dimension of the nullspace is dependent
on the number of singular values of `A` for which `iszero(a)` is true.

nu -- An Int64, which is the dimension of the nullspace.
N  -- A matrix whose columns generate the nullspace of A. Type ::Hecke.MatElem{padic}. 
"""
function nullspace(A::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    F = padic_qr(transpose(A), col_pivot=Val(true))

    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if iszero(F.R[i,:])
            push!(col_list, i)
        end
    end

    Pinv = inverse_permutation(F.p)   
    
    Q = F.Q
    inv_unit_lower_triangular!(Q)
    Qinvt = transpose(Q)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end

function nullspace(A::Hecke.SMat{padic})

    m = nrows(A)
    n = ncols(A)
    F = padic_qr(transpose(A), col_pivot=Val(true))

    # Sparse elimination sorts zero singular values to the bottom.
    # However, without column pivoting there is potential precision loss.

    
    # This function is really bad from a precision analysis standpoint. It isn't clear
    # how to deal with the insane precision gains we get from sparse elimination.
    function bad_iszero(srow)
        if iszero(srow)
            return true
        end
        return minimum(valuation.(srow.values)) >= precision(A.base_ring)
    end
    
    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if bad_iszero(F.R[i])
            push!(col_list, i)
        end
    end

    Pinv = inverse_permutation(F.p)   
    
    Q = F.Q
    println("Sparsity of Q-factor: ", sparsity(Q))
    println()
    
    badQ = matrix(Q) # The "matrix" call makes things dense. This is not ideal.
    inv_unit_lower_triangular!(badQ)   # It is probably a good idea to have a specialized QR method
                                       # that computes the inverse of Q directly, instead of Q.
    Qinvt = transpose(badQ)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end


# stable version of inverse for p-adic matrices
import Hecke.inv
"""
    inv( A::Hecke.MatElem{padic} ) -> Hecke.MatElem{padic}

Matrix inverse. Computes matrix `B` such that A*B = I, where I is the identity matrix.
Computed by solving the left-division N = M \\ I.
"""
function inv(A::Hecke.MatElem{padic})
    if size(A,1) != size(A,2)
        error("Matrix must be square.")
    end
    id = identity_matrix(A.base_ring, size(A,2))
    return rectangular_solve(A, id)
end

function inv_unit_lower_triangular!(L::Hecke.Generic.MatElem{T} where T)

    m = size(L,1)::Int64
    n = size(L,2)::Int64    
    #if !issquare(L)
    #    error("Square matrix required for inverse")
    #end
    Qp = parent(L[1,1])
    container_for_mul = Qp()
    container_for_result = Qp()
    
    for k = 1:n
        for i = k+1:n
            container_for_result=zero(Qp)
            for r=k:i-1
                Hecke.mul!(container_for_mul, L[i,r], L[r,k])
                addeq!(container_for_result,  container_for_mul)
            end
            L[i,k] = -container_for_result
        end
    end

    return
end

"""
    inv_unit_lower_triangular( A::Hecke.MatElem{padic} ) -> Hecke.MatElem{padic}

Matrix inverse, specialized to invert a lower triangular matrix with ones on the diagonal.
"""
function inv_unit_lower_triangular(L)
    L2 = deepcopy(L)
    inv_unit_lower_triangular!(L2)
    return L2
end

@doc Markdown.doc"""
    rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic}; stable::Bool=false)
                                                                        --> (nu :: Int64,N::Hecke.MatElem{padic})

Solves the linear system A*N = b. The output `nu` is the dimension of the nullspace. Parameter `stable` determines whether `padic_qr` or `svd` method is used. Default is qr (for speed).

WARNINGS:
If `A,b_input` have different precisions, maximal precision output is not guarenteed.
Underdetermined solve not implemented.
"""
function rectangular_solve(A::Hecke.MatElem{padic}, b::Hecke.MatElem{padic}; stable::Bool=false)
    return rectangular_solve!(A, deepcopy(b), stable=stable)
end

@doc Markdown.doc"""
    rectangular_solve!(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic}; stable::Bool=false)
                                                                        --> (nu :: Int64,N::Hecke.MatElem{padic})

Solves the linear system A*N = b inplace. The output `nu` is the dimension of the nullspace. Parameter `stable` determines whether `padic_qr` or `svd` method is used. Default is qr (for speed).

WARNINGS:
If `A,b_input` have different precisions, maximal precision output is not guarenteed.
Underdetermined solve not implemented.
"""
function rectangular_solve!(A::Hecke.MatElem{padic}, b::Hecke.MatElem{padic}; stable::Bool=false)
    if stable
        return _svd_rectangular_solve(A::Hecke.MatElem{padic}, b::Hecke.MatElem{padic})
    else
        return _lu_rectangular_solve(A::Hecke.MatElem{padic}, b::Hecke.MatElem{padic})
    end
end

# Specialization to lu-solve
function _lu_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        error("`A` and `b` must have the same number of rows.")
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = padic_qr(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.Q[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    println()
                    println("--- Error data: ---")
                    println("bad entry at ", i," ",j)
                    println("entries: ", b[i,j])
                    println()
                    error("Line 461: The system is inconsistent.")
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.R[i,j]
        end
        #scale = A[i, i]
        #b.row_op(i, lambda x, _: x / scale)

        if !iszero(b[i,:]) && iszero(F.R[i,i])
            println()
            println("--- Error data: ---")
            println("bad entry at row ", i)
            error("Line 480: The system is inconsistent.")
        elseif !iszero(F.R[i,i])
            b[i,:] *= inv(F.R[i,i])
        end
    end

    return b
end

# Specialization to svd-solve
function _svd_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        error("`A` and `b` must have the same number of rows.")
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = svd(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.U[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    println()
                    println("--- Error data: ---")
                    println("bad entry at ", i," ",j)
                    println("entries: ", b[i,j])
                    # println()
                    # println(b)
                    # println()
                    # println(A)
                    error("Line 533: The system is inconsistent.")
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # Scaling step
    for i in 1:n
        if !iszero(b[i,:]) && iszero(F.S[i,i])
            println()
            println("--- Error data: ---")
            error("The system is inconsistent: singular value: ", i," is zero, while `b[i,:]` is nonzero.")
        elseif !iszero(F.S[i,i])
            b[i,:] *= inv(F.S[i,i])
        end
    end

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.Vt[i,j]
        end
    end

    return b[F.q,:]
end
    
function underdetermined_solve()
    error("Not implemented.")
    return
end


#************************************************************
#
# Eigenvector iteration methods.
#
#************************************************************


#************************************************
#  Basic inverse iteration
#************************************************

# Solve for an eigenvector using inverse iteration.
# Note that the algorithm will not converge to a particular vector in general, but the norm of
#
# A*w - λ*w converges to zero. Here, λ is the unique eigenvalue closest to `shift`, (if it is unique).
#
# TODO: Separate invariant subspaces at high valuation.
const TESTFLAG=false
function inverse_iteration!(A,shift,V)

    # Note: If A is not known to precision at least one, really bad things happen.
    Qp = A.base_ring
    In = identity_matrix(A.base_ring, size(A,1))
    B  = A - shift*In
    
    if rank(B) < ncols(B)
        println("Value `shift` is exact eigenvalue. `shift` = ", shift)
        return [nullspace(B)[2]], [shift]
    end

    function normalize(V)
        maxn, m = findmax( abs.(V.entries) )
        if iszero(maxn)
            return V
        end
        return V / V[m]
    end
    
    @time pow = rectangular_solve(B,identity_matrix(B.base_ring,size(B,1)),stable=true)

    if TESTFLAG
        println("---pow---")
        println(pow)
        println("---")
        println()
    end
    
    @time for i=1:(A.base_ring.prec_max)
        V = normalize(pow*V)
        if TESTFLAG
            println(V)
            println()
        end
    end
    
    if TESTFLAG
        println("---end inv iteration---")
        println()
    end

    # Test for convergence and calculate eigenvalues,
    # if the algorithm hasn't converged, check for whether the remaining subspace can be
    # split by further iteration.
    X = try
        rectangular_solve(V, A*V, stable=true)        
    catch e
        error("Error in inverse iteration. Likely a stability issue.")
    end

    nu= trace(X)//size(X,2)
    Y = X - nu*identity_matrix(Qp,size(X,2))
    
    if iszero(Y)
        # In this case, the eigenvectors are at their maximum refinement.
        return [V],[nu]
    end

    # Since only eigenvectors (mod p) are given as the initial data, the operator Y *must* be
    # zero mod p. We scale out the denominator to try again.
    
    vals_of_Y = valuation.( Y )
    min_val = minimum(vals_of_Y)

    if min_val <=0
        error("Failure of convergence in inverse iteration.")
    end

    println("Second level iteration.")
    scale_factor = Qp(Qp.p)^Int64(-min_val)
    inv_scale_factor = Qp(Qp.p)^Int64(min_val)
    Ynew = scale_factor * Y    
    E = eigspaces(Ynew)

    return Array{typeof(V),1}([V*Esp for Esp in E.spaces]),
    Array{padic,1}([inv_scale_factor*nu for nu in E.values])
end


@doc Markdown.doc"""
    inverse_iteration(A::Hecke.MatElem{padic}, shift :: padic, v ::Hecke.MatElem{padic}  ) -> Hecke.MatElem{padic}

Iterate `v = (A-shift*I)^(-1) * v`. The inverse is cached at the beginning of the computation. The columns of the entry `v` define a subspace.

If subspace iteration does not satisfy `A*v ⊆ v`, an error is raised.
"""
function inverse_iteration(A, shift, v)
    w = deepcopy(v)
    wlist,nulist = inverse_iteration!(A,shift,w)    
    return wlist,nulist
end

@doc Markdown.doc"""
    inverse_iteration_decomposition(A::Hecke.MatElem{padic}, Amp::Hecke.MatElem{nmod_mat}  ) -> values, spaces

Return types.
    values :: Array{padic,1}
    spaces :: Array{ Hecke.MatElem{padic}, 1}

Internal function. Compute an invariant subspace decomposition of `A` using its reduction mod-p `Amp`.
Makes one call to `inverse_iteration` for each eigenvalue of `Amp`. 
"""
function inverse_iteration_decomposition(A, Amp)

    Qp = A.base_ring
    E = eigspaces(Amp)

    values_lift = fill(zero(Qp), 0)
    spaces_lift = fill(zero(parent(A)), 0)

    for i in 1:length(E.values)

        # Approximate input data
        appx_eval = Qp( lift(E.values[i]) )
        appx_espace =  matrix(Qp, lift(E.spaces[i]) )

        # Apply inverse iteration step.
        wlist,nulist = inverse_iteration(A, appx_eval, appx_espace)

        # Append refined data to the main list.
        values_lift = vcat(values_lift, nulist)
        spaces_lift = vcat(spaces_lift,  wlist)
    end

    return values_lift, spaces_lift
end


###############################################################################
#
#   Power iteration decomposition
#
###############################################################################


function power_iteration_decomposition(A, Amp)

    Qp = A.base_ring
    N = Qp.prec_max
    E = eigspaces(Amp)

    restricted_maps = Array{typeof(fill(zero(Qp), 0)),1}()
    spaces_lift = Array{typeof(fill(zero(parent(A)), 0)), 1}()

    roots_and_mults = roots_with_multiplicities(Hecke.charpoly(Amp))

    if length(E.values) > 0
        M = maximum( [ a[2] for a in roots_and_mults])
    end
    
    for i in 1:length(E.values)
       
        # Approximate input data
        appx_eval = Qp( lift(E.values[i]) )
        appx_espace =  matrix(Qp, lift(E.spaces[i]) )

        # Apply power iteration step.

        B = A - appx_eval*identity_matrix(Qp,size(A,1))

        for j=1:ceil(log2(M*N))
            B = B^2
        end

        nu,V = nullspace(B)

        if nu == 0
            display( valuation.(singular_values(A)) )
            display( valuation.(singular_values(B)) )
            error("Matrix 'B' is invertible. Iteration did not converge.")
        end
        
        X = rectangular_solve(V, A*V, stable=true)
        
        # Append refined data to the main list.
        restricted_maps = vcat(restricted_maps, [X])
        spaces_lift = vcat(spaces_lift,  [V])
    end

    return restricted_maps, spaces_lift

end

###############################################################################
#
#   "Classical Algorithm"
#
###############################################################################

function _eigenspaces_by_classical(A)
    error("Classical Algorithm not implemented in Julia. I invite the brave to implement the Zaussenhaus Round-4 algorithm.")
end

###############################################################################
#
#   Hessenberg form
#
###############################################################################

# Also return a basis by default.
@doc Markdown.doc"""
    hessenberg!(A::Hecke.Generic.Mat{T} where T <: padic; basis=Val(true)) --> nothing or B::Hecke.Generic.Mat{T}

Computes the Hessenberg form of `A` inplace. If `basis=Val(true)`, also return the matrix B such that
    AV = VB
"""
function hessenberg!(A::Hecke.Generic.Mat{T} where T <: padic; basis=Val(true))
    !issquare(A) && error("Dimensions don't match in hessenberg")
    R = base_ring(A)
    n = nrows(A)
    u = R()
    t = R()

    if basis == Val(true)
        B = identity_matrix(R, size(A,1))
    end
    
    for m = 1:n - 2

        val_list = float64_valuation.(A.entries[ m+1 : n , m ])
        minn, row_pivot_index = findmin( val_list );
        if minn==Inf continue end

        i = row_pivot_index + m;
        
        # Perform a row/column swap to move the pivot to the subdiagonal
        if i > m+1
            for j = m:n
                A[i, j], A[m+1, j] = A[m+1, j], A[i, j]
            end
            for j = 1:n
                A[j, i], A[j, m+1] = A[j, m+1], A[j, i]
            end

            if basis==Val(true)
                for j = 1:n
                    B[i, j], B[m+1, j] = B[m+1, j], B[i, j]
                end
            end
        end

        # cache the inverted pivot.
        h = -inv(A[m+1, m])

        # Perform the elimination.
        for i = m + 2:n
            if iszero(A[i, m]) continue end
            
            u = Hecke.mul!(u, A[i, m], h)

            # Row operatons
            for j = 1:n
                if j > m
                    t = Hecke.mul!(t, u, A[m+1, j])
                    A[i, j] = addeq!(A[i, j], t)
                end
                    
                if basis==Val(true)
                    t = Hecke.mul!(t, u, B[m+1,j])
                    B[i,j] = addeq!(B[i,j],t)
                end
            end
            u = -u

            # Column eliminations
            for j = 1:n
                t = Hecke.mul!(t, u, A[j, i])
                A[j, m+1] = addeq!(A[j, m+1], t)
            end
            A[i, m] = R()            
        end        
    end

    if basis==Val(true)
        return B
    else
        return
    end
end

"""
    hessenberg(A::Generic.MatrixElem{T}) where {T <: padic}
> Returns the Hessenberg form of M, i.e. an upper Hessenberg matrix
> which is similar to M. The upper Hessenberg form has nonzero entries
> above and on the diagonal and in the diagonal line immediately below the
> diagonal.
> A p-adically stable form of the algorithm is used, where pivots are 
> selected carefully.
"""
function hessenberg(A::Hecke.Generic.Mat{T} where T <: padic, basis=Val(true))
   !issquare(A) && error("Dimensions don't match in hessenberg")
   M = deepcopy(A)
   B = hessenberg!(M, basis=Val(true))
   return M, B
end


#************************************************
#  QR-iteration 
#************************************************

"""
    block_schur_form(A::Hecke.Generic.Mat{T} where T <: padic) --> B,V :: Hecke.Generic.Mat{T}

Computes the block schur form `B` of a padic matrix `A`, where the
blocks correspond to the different eigenvalues of `A modulo p`. The outputs satisfy
`AV = VB`

NOTE: 
Presently, `block_shur_form` does not attempt to further refine the blocks recursively. Theoretical
details need to be worked out to make the best practical improvements of the algorithm. 
"""
function block_schur_form(A::Hecke.Generic.Mat{T} where T <: padic)

    Qp = A.base_ring
    N = Qp.prec_max
    
    # Extract data from the reduction modulo p
    Aint  = _normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)

    B, V = hessenberg(A)
    id= identity_matrix(Qp, size(B,1))
    
    for (rt,m) in roots_with_multiplicities(chiAp)
        
        lambdaI = lift(rt)*id

        # Regarding convergence. It seems like it needs a little extra time to
        # sort the terms via permutation.
        for i in 1:N*m
            F = padic_qr(B - lambdaI)

            # Note about Julia's syntax. A[:,F.p] = A*inv(P), for a permutation P.
            B = F.R[:,F.p] * F.Q + lambdaI
            V = inv_unit_lower_triangular(F.Q)*V[F.p,:]
        end        
    end
    
    return B,V
end

@doc Markdown.doc"""
    roots_with_multiplicities(f)

Returns a list of roots of the polynomial `f`, with the multiplicities occuring in the factorization.
"""
function roots_with_multiplicities(f)
    F = Hecke.factor(f)
    return [(-g(0), m) for (g,m) in F if Hecke.degree(g) == 1]
end


function _normalize_matrix(A)

    Qp = A.base_ring
    vals_of_A = valuation.( A.entries )
    min_val = minimum(vals_of_A)

    scale_factor = Qp(Qp.p)^max(0,Int64(-min_val))
    return scale_factor * A
end

"""
    eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

Compute the eigenvectors of a padic matrix iteratively.

The `method` parameter selects the method to be used to compute the eigenvectors.
The options are:

-- "inverse"
-- "classical"
-- "power"
-- "qr"
-- "schur"

The default is `inverse`, since at the moment this is the one that is implemented.

"""

function eigspaces(A::Hecke.Generic.Mat{T} where T <: padic; method="power")

    ## Input sanitization    
    if size(A,1) != size(A,2)
        error("Input matrix must be square.")
    end

    ## Set constants
    Qp = A.base_ring

    
    if iszero(A)        
        return EigenSpaceDec(Qp, [zero(Qp)] , [identity_matrix(Qp, size(A,1))] )
    elseif  size(A,1) == 1
        return EigenSpaceDec(Qp, [A[1,1]] , [identity_matrix(Qp, size(A,1))] )
    end

    # Dispatch
    if method == "classical"
        error("Not Implemented")
        
    elseif method == "inverse"
        return  _eigenspaces_by_inverse_iteration(A)
        
    elseif method == "schur" || method == "qr"
        error("Not Implemented. However, block_shur_form is available to compute the schur form.")
        
    elseif method == "power"
        return  _eigenspaces_by_power_iteration(A)
        
    else
        error("Not Implemented")
    end

end

function _modp_charpoly_data(A::Hecke.Generic.Mat{T} where T <: padic)
    Aint  = _normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)

    return Aint, Amp, chiAp
end

function _eigenspaces_by_inverse_iteration(A::Hecke.Generic.Mat{T} where T <: padic)
    
    # Extract data from the reduction modulo p
    Qp = A.base_ring
    Aint, Amp, chiAp = _modp_charpoly_data(A)    
    factors_chiAp = Hecke.factor(chiAp)

        
    if isirreducible(chiAp)
        empty_array = Array{padic,1}()
        empty_spaces_array = Array{ Hecke.Generic.Mat{padic}, 1}()
        
        return EigenSpaceDec(Qp, empty_array , empty_spaces_array )
    end
    
    # FAILSAFE DURING DEVELOPMENT...
    # Fail automatically if there are large invariant subspaces mod p
    if any( e >= 2 for (f,e) in factors_chiAp if degree(f)==1 )
        error("Development Failsafe: Not implemented when roots are not squarefree") 
    end

    # Iteration call
    values_lift, spaces_lift = inverse_iteration_decomposition(Aint, Amp)

    return EigenSpaceDec(Qp, values_lift, spaces_lift)
    
end

function _eigenspaces_by_power_iteration(A::Hecke.Generic.Mat{T} where T <: padic)

    Qp = A.base_ring
    Aint, Amp, chiAp = _modp_charpoly_data(A)
    factors_chiAp = Hecke.factor(chiAp)
    
    empty_array = Array{padic,1}()
    empty_spaces_array = Array{ Hecke.Generic.Mat{padic}, 1}()    
        
    if isirreducible(chiAp)        
        return EigenSpaceDec(Qp, empty_array , empty_spaces_array )
    end

    factor_multiplicities = collect(Base.values(factors_chiAp.fac))

    #println(Amp)    
    #println(factor_multiplicities)
    
    # Check to ensure chiAp is not an n-th power
    if length(factors_chiAp) == 1 && factor_multiplicities[1] == size(A,1)
        return try
            _eigenspaces_by_classical(A)
        catch e
            println()
            println("WARNING: Classical Algorithm not implemented. Skipping this eigenvalue...")
            println()
            EigenSpaceDec(Qp, empty_array , empty_spaces_array )
        end
    end

    # Iteration call
    restricted_maps, invariant_blocks = power_iteration_decomposition(Aint, Amp)

    # Postprocessing
    values = fill(zero(Qp), 0)
    spaces = fill(zero(parent(Aint)), 0)    
    
    for i = 1:length(restricted_maps)

        X = restricted_maps[i]

        if size(X,1) == 1
            push!(values, X[1,1])
            push!(spaces, invariant_blocks[i])
        else
            # Recursive call
            E = _eigenspaces_by_power_iteration(X)

            #merge
            for j = 1:length(E.spaces)
                push!(values, E.values[j])
                push!(spaces, invariant_blocks[i]*E.spaces[j])
            end
            
            # for i = 1:length(E.values)
            #     push!(values, E.values[i])
            #     push!(spaces, E.spaces[i])
            # end
            
        end
        
    end
    
    return EigenSpaceDec(Qp, values, spaces)
end


############################################################################################

"""
    block_data(A)
(Non-critical testing function). Print the valuations of the main/sub diagonal.
"""
function block_data(A)

    n = size(A,2)
    data = Array{Any,2}(nothing, 2, n)
    
    for i=1:n-1
        data[:,i] = valuation.([A[i,i], A[i+1,i]])
    end

    data[:,n] = [valuation(A[n,n]), nothing]
    
    return data
end
