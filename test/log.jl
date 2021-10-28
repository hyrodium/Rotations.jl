@testset "log_3D" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 UnitQuaternion, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)
    oneaxis_types = (RotX, RotY, RotZ)

    @testset "$(T)" for T in all_types, F in (one, rand)
        R = F(T)
        @test R ≈ exp(log(R))
    end
end

@testset "log_2D" begin
    all_types = (RotMatrix{2}, Angle2d)

    @testset "$(T)" for T in all_types, θ in 0.0:0.1:π
        R = F(T)
        @test R ≈ exp(log(R))
    end
end
