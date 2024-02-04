module RotationsRecipesBaseExt

using RecipesBase
using Rotations

@recipe function f(R::Rotation{3}; origin=(0,0,0))
    e₁ = R[:,1]
    e₂ = R[:,2]
    e₃ = R[:,3]
    ox, oy, oz = origin
    @series begin
        primary := false
        color := :red
        [ox,e₁[1]+ox],[oy,e₁[2]+oy],[oz,e₁[3]+oz]
    end
    @series begin
        primary := false
        color := :green
        [ox,e₂[1]+ox],[oy,e₂[2]+oy],[oz,e₂[3]+oz]
    end
    @series begin
        primary := false
        color := :blue
        [ox,e₃[1]+ox],[oy,e₃[2]+oy],[oz,e₃[3]+oz]
    end
    seriestype := :scatter
    delete!(plotattributes, :origin)
    [ox], [oy], [oz]
end

end
