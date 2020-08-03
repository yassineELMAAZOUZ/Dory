
using Dory, Test

@testset "QR method tests" begin

    Qp = PadicField(7,20)
    A = random_test_matrix(Qp,5)

    F = padic_qr(A)
    @test A[F.p,F.q] == F.Q*F.R
    @test isupper_triangular(F.R)
    
    B = random_test_matrix(Qp, 10)
    F = padic_qr(B, col_pivot=Val(true))

    @test B[F.p,F.q] == F.Q*F.R
    @test isupper_triangular(F.R)

end


@testset "Block Schur form tests" begin
    
    Qp = PadicField(7,20)
    A = random_test_matrix(Qp,5)

    S, V = Dory.block_schur_form(A)

    @test valuation(det(V)) == 0
    @test inv(V)*S*V == A
    
    
end

nothing
