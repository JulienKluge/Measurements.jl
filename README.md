# Measurements.jl

| **Documentation**                       | [**Package Evaluator**][pkgeval-link] | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url]       | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

Introduction
------------

`Measurements.jl` is a package that allows you to define numbers with
[uncertainties](https://en.wikipedia.org/wiki/Measurement_uncertainty), perform
calculations involving them, and easily get the uncertainty of the result
according to
[linear error propagation theory](https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
This library is written in [Julia](http://julialang.org/), a modern high-level,
high-performance dynamic programming language designed for technical computing.

When used in the
[Julia interactive session](http://docs.julialang.org/en/stable/manual/getting-started/),
it can serve also as an easy-to-use calculator.

### Features List ###

* Support for most mathematical operations available in Julia standard library
  and special functions
  from [`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl)
  package, involving real and complex numbers.  All existing functions that
  accept `AbstractFloat` (and `Complex{AbstractFloat}` as well) arguments and
  internally use already supported functions can in turn perform calculations
  involving numbers with uncertainties without being redefined.  This greatly
  enhances the power of `Measurements.jl` without effort for the users
* Functional correlation between variables is correctly handled, so `x-x ≈
  zero(x)`, `x/x ≈ one(x)`, `tan(x) ≈ sin(x)/cos(x)`, `cis(x) ≈ exp(im*x)`,
  etc...
* Support for
  [arbitrary precision](http://docs.julialang.org/en/stable/manual/integers-and-floating-point-numbers/#arbitrary-precision-arithmetic)
  (also called multiple precision) numbers with uncertainties.  This is useful
  for measurements with very low relative error
* Define arrays of measurements and perform calculations with them.  Some linear
  algebra functions work out-of-the-box
* Propagate uncertainty for any function of real arguments (including functions
  based on
  [C/Fortran calls](http://docs.julialang.org/en/stable/manual/calling-c-and-fortran-code/)),
  using `@uncertain`
  [macro](http://docs.julialang.org/en/stable/manual/metaprogramming/)
* Function to get the derivative and the gradient of an expression with respect
  to one or more independent measurements
* Functions to calculate
  [standard score](https://en.wikipedia.org/wiki/Standard_score) and
  [weighted mean](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean)
* Parse strings to create measurement objects
* Easy way to attach the uncertainty to a number using the `±` sign as infix
  operator.  This syntactic sugar makes the code more readable and visually
  appealing
* Extensible in combination with external packages: you can propagate errors of
  measurements with their physical units, perform numerical integration
  with [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl), numerical and
  automatic differentiation, and much more.
* Integration with [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl).

Further features are expected to come in the future, see the section "How Can I
Help?" and the TODO list below.

The method used to handle functional correlation is described in this paper:

* M. Giordano, 2016, "Uncertainty propagation with functionally correlated
  quantities", [arXiv:1610.08716](http://arxiv.org/abs/1610.08716)
  (Bibcode:
  [`2016arXiv161008716G`](http://adsabs.harvard.edu/abs/2016arXiv161008716G))

If you use use this package for your research, please cite it.

### Documentation ###

The complete manual of `Measurements.jl` is available at
https://juliaphysics.github.io/Measurements.jl/stable/.  There, people
interested in the details of the package, in order integrate the package in
their workflow, can can find a technical appendix explaining how the package
internally works.

Installation
------------

`Measurements.jl` is available for Julia 0.7 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the commands

```julia
julia> Pkg.update()
julia> Pkg.add("Measurements")
```

Older versions of this package are also available for Julia 0.4-0.6.

Usage
-----

After installing the package, you can start using it with

```julia
using Measurements
```

The module defines a new `Measurement` data type.  `Measurement` objects can be
created with the two following constructors:

``` julia
measurement(value, uncertainty)
value ± uncertainty
```

where

* `value` is the nominal value of the measurement
* `uncertainty` is its uncertainty, assumed to be a
  [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation).

They are both subtype of `AbstractFloat`.  Some keyboard layouts provide an easy
way to type the `±` sign, if your does not, remember you can insert it in Julia
REPL with `\pm` followed by `TAB` key.  You can provide `value` and
`uncertainty` of any subtype of `Real` that can be converted to `AbstractFloat`.
Thus, `measurement(42, 33//12)` and `pi ± 0.1` are valid.

`measurement(value)` creates a `Measurement` object with zero uncertainty, like
mathematical constants.  See below for further examples.

Every time you use one of the constructors above, you define a *new independent*
measurement.  Instead, when you perform mathematical operations involving
`Measurement` objects you create a quantity that is not independent, but rather
depends on really independent measurements.

Most mathematical operations are instructed, by
[operator overloading](https://en.wikipedia.org/wiki/Operator_overloading), to
accept `Measurement` type, and uncertainty is calculated exactly using analityc
expressions of functions’ derivatives.

In addition, it is possible to create a `Complex` measurement with
`complex(measurement(a, b), measurement(c, d))`.

``` julia
measurement(string)
```

`measurement` function has also a method that enables you to create a
`Measurement` object from a string.

This module extends many methods defined in Julia’s mathematical standard
library, and some methods from widespread third-party packages as well.  This is
the case for most special functions
in [`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl)
package, and the `quadgk` integration routine
from [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) package.  See the
full manual for details.

Examples
--------

``` julia
julia> using Measurements

julia> a = measurement(4.5, 0.1)
4.5 ± 0.1

julia> b = 3.8 ± 0.4
3.8 ± 0.4

julia> 2a + b
12.8 ± 0.4472135954999579

julia> x = 8.4 ± 0.7

julia> x - x
0.0 ± 0.0

julia> x/x
1.0 ± 0.0

julia> x*x*x - x^3
0.0 ± 0.0

julia> sin(x)/cos(x) - tan(x)
-2.220446049250313e-16 ± 0.0 # They are equal within numerical accuracy
```

License
-------

The `Measurements.jl` package is licensed under the MIT "Expat" License.  The
original author is Mosè Giordano.

Please, cite the paper Giordano 2016 (http://arxiv.org/abs/1610.08716) if you
employ this package in your research work.


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://juliaphysics.github.io/Measurements.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliaphysics.github.io/Measurements.jl/stable/

[pkgeval-link]: http://pkg.julialang.org/?pkg=Measurements

[pkg-0.5-img]: http://pkg.julialang.org/badges/Measurements_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/detail/Measurements.html
[pkg-0.6-img]: http://pkg.julialang.org/badges/Measurements_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/detail/Measurements.html

[travis-img]: https://travis-ci.org/JuliaPhysics/Measurements.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaPhysics/Measurements.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/u8mg5dlhyb1vjcpe?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/measurements-jl

[coveral-img]: https://coveralls.io/repos/github/JuliaPhysics/Measurements.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/JuliaPhysics/Measurements.jl?branch=master

[codecov-img]: https://codecov.io/gh/JuliaPhysics/Measurements.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaPhysics/Measurements.jl
