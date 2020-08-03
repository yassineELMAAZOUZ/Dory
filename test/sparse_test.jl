
# Some sparse linear algebra tests.

using Dory

# generate sparse matrix over Q
p = 11
Qp = PadicField(p,3)


function random_sparse_matrix(Qp)
    A = sparse_matrix(Qp)

    # adds sparse rows with three entries each
    for i=1:10
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

for j=1:100

    println("\n\n============================================")
    println( "NEW TEST ")
    println("===\n")
    
    A = random_sparse_matrix(Qp)

    display( Dory.float64_valuation.(Matrix(A) ) )
    
    @time data = padic_qr(A)

    display( Dory.float64_valuation.(Matrix(data.R) ) )
end


