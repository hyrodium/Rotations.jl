@testset "log" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 QuatRotation, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY,
                 RotMatrix{2}, Angle2d)

    @testset "$(T)" for T in all_types, F in (one, rand)
        R = F(T)
        @test R ≈ exp(log(R))
        @test log(R) isa InfinitesimalRotation
        @test exp(log(R)) isa Rotation
    end
end

@testset "exp(zero)" begin
    all_types = (InfinitesimalRotMatrix{3}, InfinitesimalRotationVec,
                 InfinitesimalRotMatrix{2}, InfinitesimalAngle2d)

    @testset "$(T)" for T in all_types
        r = zero(T)
        @test one(exp(r)) ≈ exp(r)
        @test exp(r) isa Rotation
    end
end
