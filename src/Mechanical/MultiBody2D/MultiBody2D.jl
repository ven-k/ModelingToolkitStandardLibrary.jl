module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ...DynamicQuantities: @u_str
using ..TranslationalPosition

@parameters t [unit = u"s"]
D = Differential(t)

export Link
include("components.jl")

end
