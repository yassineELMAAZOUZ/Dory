
using Dory, Random

# Enforce that tests are identical.
Random.seed!(123)
Qp = PadicField(41,100)

function rand_input(n)
    return matrix(Qp, n, n, [rand_padic_int(Qp) for i=1:n for j=1:n])
end

# Run the first time to compile
A = rand_input(5)
E  = eigspaces(A)
E2 = Dory.block_schur_form(A)


# timings
num_samp = 10
for n in [10,100,200]

    tot_time_power = 0
    tot_time_class = 0
    tot_time_schur = 0
    
    for j=1:num_samp
        A = rand_input(n)
    
        t0 = time()
        E = eigspaces(A)
        t1 = time()
        tot_time_power += t1-t0
        
        t0 = time()
        E2 = Dory.block_schur_form(A)
        t1 = time()
        tot_time_schur += t1-t0
        
    end

    avg_time_power = tot_time_power/num_samp
    avg_time_schur = tot_time_schur/num_samp

    
    @info "" n avg_time_power avg_time_schur
    display([avg_time_power, avg_time_schur])
    
end
