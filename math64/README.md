# math

[![Go Reference](https://pkg.go.dev/badge/github.com/vkngwrapper/math.svg)](https://pkg.go.dev/github.com/vkngwrapper/math)

`go get github.com/vkngwrapper/math`

`vkngwrapper/math` is a fast, low-allocation 3d math library, oriented toward Vulkan. That means that in 
 projection-focused methods, it uses a right-handed coordinate space and a clip space with a 0-1 Z-axis. It
 is developed to be used with `vkngwrapper`, but can be used for any vulkan or vulkan-like 3d space which
 requires matrix or quaternion arithmetic.

## How Fast Is It?

I ran four benchmarks comparing vkng math to three other pure go 3d math libraries: 
 [g3n math32](https://github.com/g3n/engine/tree/master/math32), 
 [go3d](https://github.com/ungerik/go3d), and [go-gl/mathgl](https://github.com/go-gl/mathgl).
 The four benchmarks were the following:

* Build a transform matrix that performs a translate, a rotate, then a scale, and then transform a 
  vector through the matrix. Do this in the fastest way available in the library.
* Rotate a vector around an axis
* Rotate a vector with a quaternion
* Build matrices that individually perform a translate, a rotate, and a scale, multiply them together
  into a single matrix, and then transform a vector through the matrix.

Results:

```shell
$ go test -test.bench=. .
goos: windows
goarch: amd64
pkg: github.com/vkngwrapper/math
cpu: AMD Ryzen 7 5800X 8-Core Processor
BenchmarkTransformVec4G3n-16                    23740468                50.53 ns/op            0 B/op          0 allocs/op
BenchmarkTransformVec4MathGL-16                 10808385               112.9 ns/op             0 B/op          0 allocs/op
BenchmarkTransformVec4Go3D-16                   28565443                40.41 ns/op            0 B/op          0 allocs/op
BenchmarkTransformVec4VkngMath-16               34278466                34.90 ns/op            0 B/op          0 allocs/op
BenchmarkRotateVec3G3n-16                       49986253                23.85 ns/op            0 B/op          0 allocs/op
BenchmarkRotateVec3MathGL-16                    22217284                54.02 ns/op            0 B/op          0 allocs/op
BenchmarkRotateVec3Go3D-16                      42848725                28.71 ns/op            0 B/op          0 allocs/op
BenchmarkRotateVec3VkngMath-16                  59987402                20.68 ns/op            0 B/op          0 allocs/op
BenchmarkRotateQuaternionG3n-16                 24484348                49.14 ns/op            0 B/op          0 allocs/op
BenchmarkRotateQuaternionMathGL-16              19043356                62.61 ns/op            0 B/op          0 allocs/op
BenchmarkRotateQuaternionGo3D-16                11108570               105.7 ns/op             0 B/op          0 allocs/op
BenchmarkRotateQuaternionVkngMath-16            25526482                46.75 ns/op            0 B/op          0 allocs/op
BenchmarkMatrixMultTransformG3n-16              20334501                59.86 ns/op            0 B/op          0 allocs/op
BenchmarkMatrixMultTransformMatGL-16            16434730                72.91 ns/op            0 B/op          0 allocs/op
BenchmarkMatrixMultTransformGo3D-16             15996331                74.72 ns/op            0 B/op          0 allocs/op
BenchmarkMatrixMultTransformVkngMath-16         21813343                54.29 ns/op            0 B/op          0 allocs/op
PASS
ok      github.com/vkngwrapper/math     20.310s
```

![the results in graph form](https://github.com/vkngwrapper/math/blob/main/readme/chart.png?raw=true)

As you can see, vkng math is faster than the alternatives, sometimes by quite a bit. If you are looking for a
 fast 3d math library that is more opengl friendly, I would recommend g3n's math library given these results.

## Why not SIMD?

Originally, this library had been planned to include SIMD support with the intent of adding NEON later on.
 However, when I began to integrate SIMD into the library, it became clear that Go has a (small, single-digit
 nanoseconds) overhead when calling into asm methods that for the most part cancelled out gains from SIMD,
 even for relatively hefty methods like matrix multiply. As a result, I chose to leave out SIMD support
 unless the situation changes.

## How to use it?

This library is relatively straightforward- vkng math objects are intended to be used on the stack where
 possible, in order to avoid allocations. One might build a relatively complicated matrix like this:

```go
		var mat Mat4x4[float32]
		mat.SetTranslation(1, 1, 1)
		mat.RotateY(1)
		mat.Scale(1.5, 1.5, 1.5)
```

Working with a new matrix or quaternion should always start with one of the `Set*` methods (`SetIdentity`, at least)
 and further modifications can be made to it with one of the non-`Set*` methods. One might use the above matrix 
 like this:

```go
		v := Vec4[float32]{5, 10, 15, 1}
		v.Transform(&mat)
```

The above code makes zero allocations and can complete in approximately 35ns. See the godoc reference linked
 at the top of this README to learn about all the object types and methods available in this library.
