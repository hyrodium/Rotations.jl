# 2D Rotation Generator

```@setup generator2d
using Rotations
```

## `RotMatrixGenerator2`

**example**

```@repl generator2d
m = rand(2,2)
s = RotMatrixGenerator{2}(m - m')
exp(s)
log(exp(s))
```

## `Angle2dGenerator`

**example**

```@repl generator2d
s = Angle2dGenerator(0.42)
exp(s)
log(exp(s))
rotation_angle(exp(s))
rotation_angle(exp(2s))
rotation_angle(exp(2s)) / rotation_angle(exp(s))
```
