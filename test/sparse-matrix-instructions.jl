
# generate sparse matrix over Q
S = sparse_matrix(QQ)

# adds two sparse rows with two entries
push!(S, sparse_row(QQ, [1,100], [Q(1), Q(2)]))
push!(S, sparse_row(QQ, [2,100], [Q(17), Q(2)]))

# read data from the sparse matrix
for R in S.rows
    for (p,v) = R
        println("in pos ", p, " we have ", v)
    end
end

# basic gauss step, adds 4*S[1] to S[2] in place (changes S)
Hecke.add_scaled_row!(S, 1, 2, QQ(4))
