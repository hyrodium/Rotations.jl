module RotationsRecipesBaseExt

using RecipesBase
using Rotations

@recipe function f(R::Rotation{3}; origin=SVector(0,0,0), boxsize=0.2)
    r = boxsize
    e₁ = R[:,1]
    e₂ = R[:,2]
    e₃ = R[:,3]
    ox, oy, oz = origin
    ps = vec([SVector(ox,oy,oz)+R*SVector(x,y,z) for x in (-r,r), y in (-r,r), z in (-r,r)])
    xs = getindex.(ps,1)
    ys = getindex.(ps,2)
    zs = getindex.(ps,3)
    @series begin
        primary := false
        color := RGB(1,0,0)
        [e₁[1]*r+ox,e₁[1]+ox],[e₁[2]*r+oy,e₁[2]+oy],[e₁[3]*r+oz,e₁[3]+oz]
    end
    @series begin
        primary := false
        color := RGB(0,1,0)
        [e₂[1]*r+ox,e₂[1]+ox],[e₂[2]*r+oy,e₂[2]+oy],[e₂[3]*r+oz,e₂[3]+oz]
    end
    @series begin
        primary := false
        color := RGB(0,0,1)
        [e₃[1]*r+ox,e₃[1]+ox],[e₃[2]*r+oy,e₃[2]+oy],[e₃[3]*r+oz,e₃[3]+oz]
    end
    seriestype := :mesh3d
    connections := (
        # Somehow 0-based indexing
        # https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref047/
        [1-1,1-1,1-1,1-1,1-1,1-1,8-1,8-1,8-1,8-1,8-1,8-1,],
        [2-1,6-1,4-1,3-1,7-1,5-1,4-1,3-1,7-1,5-1,6-1,2-1,],
        [6-1,5-1,2-1,4-1,3-1,7-1,3-1,7-1,5-1,6-1,2-1,4-1,],
    )
    # This connections can be 1-based indexing, but this throws an error on PythonPlot.
    # https://discourse.julialang.org/t/how-to-plot-a-cube-in-3d-in-plots-jl/86919/2?u=hyrodium
    # connections := [(1,2,6),(1,6,5),(1,4,2),(1,3,4),(1,7,3),(1,5,7),(8,4,3),(8,3,7),(8,7,5),(8,5,6),(8,6,2),(8,2,4)]
    delete!(plotattributes, :origin)
    delete!(plotattributes, :boxsize)
    xs,ys,zs
end

end
