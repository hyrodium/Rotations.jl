@testset "nearest_rotation" begin
    @testset "$(N)-dim" for N in [2,3]
        for _ in 1:100
            M = randn(N,N)
            R = nearest_rotation(M)
            V = M/R
            U = R\M

            # R is rotation matrix
            @test isrotation(R)
            # U, V are (approximately) symmetric matrices
            @test U ≈ U'
            @test V ≈ V'
            # Eigen values are non-negative, except for the first value.
            @test all(eigvals(U)[2:end] .≥ 0)
            @test all(eigvals(V)[2:end] .≥ 0)
        end

        for _ in 1:100
            M = @SMatrix randn(N,N)
            R = nearest_rotation(M)
            V = M/R
            U = R\M

            @test isrotation(R)
            @test U ≈ U'
            @test V ≈ V'
            # Can't calculate eigen values of general SMatrix.
            @test all(eigvals(Symmetric(U))[2:end] .≥ 0)
            @test all(eigvals(Symmetric(V))[2:end] .≥ 0)
        end
    end
end
