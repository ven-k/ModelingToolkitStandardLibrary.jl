"""
    Gain(k; name)

Output the product of a gain value with the input signal.

# Parameters:

  - `k`: Scalar gain

# Connectors:

  - `input`
  - `output`
"""
@model Gain begin # (k; name)
    @extend u, y = siso = SISO()
    @parameters begin
        k=k, [description = "Gain of Gain"]
    end
    @equations begin
        y ~ k * u
    end
end

"""
    MatrixGain(K::AbstractArray; name)

Output the product of a gain matrix with the input signal vector.

# Parameters:

  - `K`: Matrix gain

# Connectors:

  - `input`
  - `output`
"""
@model MatrixGain  begin#(K::AbstractArray; name)
    begin
        nout, nin = size(K, 1), size(K, 2)
    end
    @extend begin
        input = RealInput(; nin = nin)
        output = RealOutput(; nout = nout)
    end
    @equations begin
        [output.u[i] ~ sum(K[i, j] * input.u[j] for j in 1:nin) for i in 1:nout]
    end
end

"""
    Sum(n::Int; name)

Output the sum of the elements of the input port vector.

# Parameters:

  - `n`: Input port dimension

# Connectors:

  - `input`
  - `output`
"""
## TODO: compose
@model Sum begin#(n::Int; name)
    @extend begin
        input = RealInput(; nin = n)
        output = RealOutput()
    end
    @equations begin
        output.u ~ sum(input.u)
    end
    compose(ODESystem(eqs, t, [], []; name = name), [input, output])
end

"""
    Feedback(;name)

Output difference between reference input (input1) and feedback input (input2).

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@model Feedback begin
    @extend begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end
    @equations begin
        output.u ~ input1.u - input2.u
    end
end

"""
    Add(;name, k1=1, k2=1)

Output the sum of the two scalar inputs.

# Parameters:

  - `k1`: Gain for first input
  - `k2`: Gain for second input

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@model Add begin # (; name, k1 = 1, k2 = 1)
    @extend begin
        input1 = RealInput()
        input2 = RealInput()
        output = RealOutput()
    end
    @parameters begin
        k1=k1, [description = "Gain of Add $name input1"]
        k2=k2, [description = "Gain of Add $name input2"]
    end
    @equations begin
        output.u ~ k1 * input1.u + k2 * input2.u
    end
end

"""
    Add(;name, k1=1, k2=1,k3=1)

Output the sum of the three scalar inputs.

# Parameters:

  - `k1`: Gain for first input
  - `k2`: Gain for second input
  - `k3`: Gain for third input

# Connectors:

  - `input1`
  - `input2`
  - `input3`
  - `output`
"""
@model Add3(; name, k1 = 1, k2 = 1, k3 = 1)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named input3 = RealInput()
    @named output = RealOutput()
    pars = @parameters(k1=k1, [description = "Gain of Add $name input1"],
        k2=k2, [description = "Gain of Add $name input2"],
        k3=k3, [description = "Gain of Add $name input3"],)
    eqs = [
        output.u ~ k1 * input1.u + k2 * input2.u + k3 * input3.u,
    ]
    return compose(ODESystem(eqs, t, [], pars; name = name), input1, input2, input3, output)
end

"""
    Product(;name)

Output product of the two inputs.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@model Product(; name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs = [
        output.u ~ input1.u * input2.u,
    ]
    return compose(ODESystem(eqs, t, [], []; name = name), input1, input2, output)
end

"""
    Division(;name)

Output first input divided by second input.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@model Division(; name)
    @named input1 = RealInput()
    @named input2 = RealInput(u_start = 1.0) # denominator can not be zero
    @named output = RealOutput()
    eqs = [
        output.u ~ input1.u / input2.u,
    ]
    return compose(ODESystem(eqs, t, [], []; name = name), input1, input2, output)
end

"""
    StaticNonLinearity(func ;name)

Applies the given function to the input.

If the given function is not composed of simple core methods (e.g. sin, abs, ...), it has to be registered via `@register_symbolic func(u)`

# Connectors:

  - `input`
  - `output`
"""
@model StaticNonLinearity(func; name)
    @named siso = SISO()
    @unpack u, y = siso
    eqs = [y ~ func(u)]
    extend(ODESystem(eqs, t, [], []; name = name), siso)
end

"""
    Abs(;name)

Output the absolute value of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Abs(; name) = StaticNonLinearity(abs; name)

"""
    Sign(;name)

Output the sign of the input

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sign(; name) = StaticNonLinearity(sign; name)

"""
    Sqrt(;name)

Output the square root of the input (input >= 0 required).

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sqrt(; name) = StaticNonLinearity(sqrt; name)

"""
    Sin(;name)

Output the sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sin(; name) = StaticNonLinearity(sin; name)

"""
    Cos(;name)

Output the cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Cos(; name) = StaticNonLinearity(cos; name)

"""
    Tan(;name)

Output the tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Tan(; name) = StaticNonLinearity(tan; name)

"""
    Asin(;name)

Output the arc sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Asin(; name) = StaticNonLinearity(asin; name)

"""
    Acos(;name)

Output the arc cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Acos(; name) = StaticNonLinearity(acos; name)

"""
    Atan(;name)

Output the arc tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Atan(; name) = StaticNonLinearity(atan; name)

"""
    Atan2(;name)

Output the arc tangent of the input.

# Connectors:

  - `input1`
  - `input2`
  - `output`
"""
@model Atan2(; name)
    @named input1 = RealInput()
    @named input2 = RealInput()
    @named output = RealOutput()
    eqs = [
        output.u ~ atan(input1.u, input2.u),
    ]
    compose(ODESystem(eqs, t, [], []; name = name), [input1, input2, output])
end

"""
    Sinh(;name)

Output the hyperbolic sine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Sinh(; name) = StaticNonLinearity(sinh; name)

"""
    Cosh(;name)

Output the hyperbolic cosine of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Cosh(; name) = StaticNonLinearity(cosh; name)

"""
    Tanh(;name)

Output the hyperbolic tangent of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Tanh(; name) = StaticNonLinearity(tanh; name)

"""
    Exp(;name)

Output the exponential (base e) of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Exp(; name) = StaticNonLinearity(exp; name)

"""
    Log(;name)

Output the natural (base e) logarithm of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Log(; name) = StaticNonLinearity(log; name)

"""
    Log10(;name)

Output the base 10 logarithm of the input.

# Connectors:

See [`StaticNonLinearity`](@ref)
"""
@component Log10(; name) = StaticNonLinearity(log10; name)
