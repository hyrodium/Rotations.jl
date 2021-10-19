using Rotations
using Documenter

DocMeta.setdocmeta!(Rotations, :DocTestSetup, :(using Rotations); recursive=true)

makedocs(;
    modules=[Rotations],
    repo="https://github.com/JuliaGeometry/Rotations.jl/blob/{commit}{path}#{line}",
    sitename="Rotations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeometry.github.io/Rotations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Types" => "types.md",
        "2D rotations" => [
            "Matrix" => "2d_matrix.md",
            "Angle" => "2d_angle.md",
            ],
        "3D rotations" => [
            "Matrix" => "3d_matrix.md",
            "Euler angles" => "3d_euler.md",
            "Angle and axis" => "3d_angleaxis.md",
            "Quaternion" => "3d_quaternion.md",
            "MRP" => "3d_mrp.md",
            ],
        "Useful functions" => "functions.md",
        "Random" => "random.md",
        "Function Reference" => "functionreference.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeometry/Rotations.jl",
)
