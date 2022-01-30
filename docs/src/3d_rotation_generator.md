# 3D Rotation Generator

```@setup generator3d
using Rotations
```

## `RotMatrixGenerator3`

**example**

```@repl generator3d
m = rand(3,3)
s = RotMatrixGenerator{3}(m - m')
exp(s)
```

## `RotationVecGenerator`

**example**

```@repl generator3d
s = RotationVecGenerator(0.1,0.2,0.3)
exp(s)
log(exp(s))
rotation_angle(exp(s))
rotation_angle(exp(2s))
rotation_angle(exp(2s)) / rotation_angle(exp(s))
```
