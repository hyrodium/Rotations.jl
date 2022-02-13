# Rotation Generator Types

## Abstract rotation generators
A matrix ``R`` is called *skew-symmetric matrix* if ``R`` satisfies

```math
\begin{aligned}
R^\top &= -R.
\end{aligned}
```

In `Rotations.jl`, there's an abstract type for skew-symmetric matrix, `RotationGenerator{L}`.
Where `L` is a size of the skew-symmetric matrix.

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
    * Skew symmetric matrix in 2 dimensional Euclidean space.
* `Angle2dGenerator`
    * Parametrized with one real number like `Angle2d`.

### 3D rotations
* `RotMatrixGenerator3{T}`
    * Skew symmetric matrix in 3 dimensional Euclidean space.
* `RotationVecGenerator`
    * Rotation generator around given axis. The length of axis vector represents its angle.
