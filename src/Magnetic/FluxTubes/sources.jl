"""
Constant magnetomotive force.

Parameters:

  - `V_m`: [A] Magnetic potential difference
"""
@mtkmodel ConstantMagneticPotentialDifference begin
  @extend Phi, V_m = twoport = TwoPort(; V_m = 0.0)
end
ConstantMagneticPotentialDifference.f(; V_m, name) = ConstantMagneticPotentialDifference.f(; twoport__V_m = V_m, name)

"""
Source of constant magnetic flux.

Parameters:

  - `Phi`: [Wb] Magnetic flux
"""
@mtkmodel ConstantMagneticFlux begin
  @extend Phi, V_m = twoport = TwoPort(; Phi = 1.0)
end
ConstantMagneticFlux.f(; Phi, name) = ConstantMagneticFlux.f(; twoport__Phi = Phi, name)
