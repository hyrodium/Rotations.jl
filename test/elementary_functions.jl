@testset "log" begin
    all_types = (RotMatrix3, RotMatrix{3}, AngleAxis, RotationVec,
                 QuatRotation, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY,
                 RotMatrix2, RotMatrix{2}, Angle2d)

    @testset "$(T)" for T in all_types, F in (one, rand)
        R = F(T)
        @test R ≈ exp(log(R))
        @test log(R) isa RotationGenerator
        @test exp(log(R)) isa Rotation
    end

    @testset "$(N)-dim" for N in 1:5
        M = @SMatrix rand(N,N)
        R = nearest_rotation(M)
        @test isrotationgenerator(log(R))
        @test log(R) isa RotMatrixGenerator
        @test exp(log(R)) isa RotMatrix
    end
end

@testset "exp(zero)" begin
    all_types = (RotMatrixGenerator{3}, RotationVecGenerator,
                 RotMatrixGenerator{2}, Angle2dGenerator)

    @testset "$(T)" for T in all_types
        r = zero(T)
        @test one(exp(r)) ≈ exp(r)
        @test exp(r) isa Rotation
    end
end

@testset "exp(::RotMatrixGenerator)" begin
    for N in 2:3
        r = zero(RotMatrixGenerator{N})
        @test r isa RotMatrixGenerator{N}
        @test exp(r) isa RotMatrix{N}
    end
end
