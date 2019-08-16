
# Simple tests for now.

using Dory, Test

Qp = PadicField(7,20)
A = random_test_matrix(Qp,5)

F = padic_qr(A)
@test A[F.p,F.q] == F.Q*F.R

B = random_test_matrix(Qp, 10)
F = padic_qr(B, col_pivot=Val(true))

@test B[F.p,F.q] == F.Q*F.R
