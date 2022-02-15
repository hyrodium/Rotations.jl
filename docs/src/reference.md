# Reference for rotations

## Published articles
* ["Fundamentals of Spacecraft Attitude Determination and Control" by Markley and Crassidis](https://link.springer.com/book/10.1007/978-1-4939-0802-8)
    * See sections 2.1-2.9 for 3d rotation parametrization.
    * See sections 3.1-3.2 for kinematics.
* ["Derivation of the Euler-Rodrigues formula for three-dimensional rotations from the general formula for four-dimensional rotations" by Johan Ernest Mebius](https://arxiv.org/abs/math/0701759)
    * The conversion `RotMatrix{3}` to `QuatRotaiton` is based on this paper.
* ["Computing Exponentials of Skew Symmetric Matrices And Logarithms of Orthogonal Matrices" by Jean Gallier and Dianna Xu](https://cs.brynmawr.edu/~dxu/206-2550-2.pdf)
    * See this article for log and exp of rotation matrices in higher dimensions.

## Wikipedia articles
* [Rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix)
* [Quaternion](https://en.wikipedia.org/wiki/Quaternion)
* [Wahba's problem](https://en.wikipedia.org/wiki/Wahba%27s_problem)
    * Wahba's problem is related to the `nearest_rotation` function.
* [Polar decomposition](https://en.wikipedia.org/wiki/Polar_decomposition)
    * Polar decomposition is also related to the `nearest_rotation` function.

## Others
* ["Quaternions and 3d rotation, explained interactively" by 3Blue1Brown](https://www.youtube.com/watch?v=zjMuIxRvygQ)
* ["Visualizing quaternions (4d numbers) with stereographic projection" by 3Blue1Brown](https://www.youtube.com/watch?v=d4EgbgTm0Bg)
