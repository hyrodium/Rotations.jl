@testset "plot" begin
    dir_out = joinpath(@__DIR__, "out_plot")
    rm(dir_out; force=true, recursive=true)
    mkpath(dir_out)

    @testset "2d" begin
        pl = plot(Angle2d(0.2); aspectratio=1)
        plot!(pl, Angle2d(0.3); aspectratio=1)

        path_img = joinpath(dir_out, "2d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end

    @testset "3d" begin
        N = 12
        Random.seed!(42)
        R1 = rand(QuatRotation)
        R2 = rand(QuatRotation)
        pl = plot(R1; linewidth=3, label="R1")
        plot!(pl, R2; linewidth=3, label="R2", origin=(1,0,0))
        for i in 2:N-1
            t = i/N
            plot!(pl, slerp(R1,R2,t); linewidth=3, label=nothing, origin=(t,0,0), boxsize=0.0, axissize=3/4)
        end

        path_img = joinpath(dir_out, "3d.png")
        @test !isfile(path_img)
        savefig(pl, path_img)
        @test isfile(path_img)
    end
end
