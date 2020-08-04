
using Dory, Test

@testset "QR method tests" begin

    Qp = PadicField(7,20)

    @testset "Basic version" begin
        for n=1:10
            A = zero_matrix(Qp,n)
            F = padic_qr(A)
            @test A[F.p,F.q] == F.Q*F.R
            @test isupper_triangular(F.R)

            for i=1:5
                B = random_test_matrix(Qp, n)
                F = padic_qr(B)
                @test B[F.p,F.q] == F.Q*F.R
                @test isupper_triangular(F.R)
            end
        end
    end

    @testset "Column pivoting" begin
        for n=1:10
            A = zero_matrix(Qp,n)
            F = padic_qr(A, col_pivot=Val(true))
            @test A[F.p,F.q] == F.Q*F.R
            @test isupper_triangular(F.R)

            for i=1:5
                B = random_test_matrix(Qp, n)
                F = padic_qr(B, col_pivot=Val(true))
                @test B[F.p,F.q] == F.Q*F.R
                @test isupper_triangular(F.R)
            end
        end
    end
    
end


@testset "Block Schur form tests" begin
    
    Qp = PadicField(7,30)

    for n=1:10
        for i=1:5
            A = random_test_matrix(Qp,n)
            S, V = Dory.block_schur_form(A)

            @test valuation(det(V)) == 0
            @test inv(V)*S*V == A
            @test Dory.isweak_block_schur_hessenberg(S)

            # if !(Dory.isweak_block_schur_hessenberg(S))
            #     println()
            #     display(factor(charpoly(S)))
            #     display(factor(charpoly(modp.(S))))
            #     display(valuation.(S))
            # end
        end
    end
end

nothing
