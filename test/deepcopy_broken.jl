
using Hecke

# generate sparse matrix over Q
p = 11
Qp = PadicField(p,3)

function rand_padic_int(Qp)
    return Qp(rand(1:10))
end

function random_sparse_matrix(Qp)
    A = sparse_matrix(Qp)

    # adds sparse rows with three entries each
    for i=1:7
        a,b,c = 0,0,0
        while !( a < b < c)
            a,b,c = rand(1:10), rand(1:10), rand(1:10)
        end
        
        poses = [a,b,c]
        entries= [rand_padic_int(Qp) for j=1:3]
        push!(A, sparse_row(Qp, poses, entries))
    end
    return A
end

A = random_sparse_matrix(Qp)

# This is a little inconvenient.
B = deepcopy(A)
display(A == B)
display(A.base_ring == B.base_ring)
