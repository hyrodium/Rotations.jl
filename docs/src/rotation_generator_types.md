# Rotation Generator Types

## Abstract rotation generators
A matrix ``R`` is called *rotation generator matrix* if ``R`` satisfies

```math
\begin{aligned}
R^\top &= -R^{-1}.
\end{aligned}
```

In `Rotations.jl`, there's an abstract type for rotations generator matrix, `RotationGenerator{L}`.
Where `L` is a size of the rotation generator matrix.

## Type hierarchy

```@setup hierarchy
using InteractiveUtils
```

```@repl hierarchy
using Rotations, StaticArrays
RotationGenerator <: StaticMatrix <: AbstractMatrix
subtypes(RotationGenerator{2})
subtypes(RotationGenerator{3})
```

## Overview of each type
For more information, see the sidebar page.

### 2D rotations
* `RotMatrixGenerator2{T}`
    * Rotation matrix in 2 dimensional Euclidean space.
* `Angle2dGenerator`
    * Parametrized with rotational angle.

### 3D rotations
* `RotMatrixGenerator3{T}`
    * Rotation matrix in 3 dimensional Euclidean space.
* `RotationVecGenerator`
    * Rotation around given axis. The length of axis vector represents its angle.
