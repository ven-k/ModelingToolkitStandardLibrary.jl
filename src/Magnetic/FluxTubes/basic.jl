"""
    Ground(;name)

Zero magnetic potential.
"""
@mtkmodel Ground begin
    @components begin
        port = PositiveMagneticPort()
    end
    @equations begin
        port.V_m ~ 0
    end
end

"""
    Idle(;name)

Idle running branch.
"""
@mtkmodel Idle begin
    @extend (Phi,) = two_port = TwoPort()
    @equations begin
        Phi ~ 0
    end
end

"""
    Short(;name)

Short cut branch.
"""
@mtkmodel Short begin
    @extend (V_m, ) = two_port = TwoPort()
    @equations begin
        V_m ~ 0
    end
end

"""
    Crossing(;name)

Crossing of two branches.

This is a simple crossing of two branches. The ports port_p1 and port_p2 are connected, as well as port_n1 and port_n2.
"""
@mtkmodel Crossing begin
    @components begin
        port_p1 = PositiveMagneticPort()
        port_p2 = PositiveMagneticPort()
        port_n1 = NegativeMagneticPort()
        port_n2 = NegativeMagneticPort()
    end
    @equations begin
        connect(port_p1, port_p2)
        connect(port_n1, port_n2)
    end
end

"""
    ConstantPermeance(;name, G_m=1.0)

Constant permeance.

# Parameters:

  - `G_m`: [H] Magnetic permeance
"""
@mtkmodel ConstantPermeance begin
    @extend V_m, Phi = two_port = TwoPort()
    @parameters begin
        G_m = 1.0
    end
    @equations begin
        Phi ~ G_m * V_m
    end
end

"""
    ConstantReluctance(;name, R_m=1.0)

Constant reluctance.

# Parameters:

  - `R_m`: [H^-1] Magnetic reluctance
"""
@mtkmodel ConstantReluctance begin
    @extend V_m, Phi = two_port = TwoPort(; Phi = 0.0)
    @parameters begin
        R_m = 1.0
    end
    @equations begin
        V_m ~ Phi * R_m
    end
end

"""
    ElectroMagneticConverter(;name, N, Phi = 0.0)

Ideal electromagnetic energy conversion.

The electromagnetic energy conversion is given by Ampere's law and Faraday's law respectively
V_m = N * i
N * dΦ/dt = -v

# Parameters:

  - `N`: Number of turns
  - `Phi_start`: [Wb] Initial magnetic flux flowing into the port_p
"""
@mtkmodel ElectroMagneticConverter begin
    @components begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
        p = Pin()
        n = Pin()
    end

    @variables begin
        v(t)
        i(t)
        V_m(t)
        Phi(t) = 0.0
    end
    @parameters begin
        N
    end
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        V_m ~ port_p.V_m - port_n.V_m
        0 ~ port_p.Phi + port_n.Phi
        Phi ~ port_p.Phi
    #converter equations:
        V_m ~ i * N # Ampere's law
        D(Phi) ~ -v / N
    end
end

"""
    EddyCurrent(;name, rho = 0.098e-6, l = 1, A = 1, Phi = 0.0)

For modelling of eddy current in a conductive magnetic flux tube.

# Parameters:

  - `rho`: [ohm * m] Resistivity of flux tube material (default: Iron at 20degC)
  - `l`: [m] Average length of eddy current path
  - `A`: [m^2] Cross sectional area of eddy current path
  - `Phi`: [Wb] Initial magnetic flux flowing into the port_p
"""
@mtkmodel EddyCurrent begin
    @extend (V_m, Phi) = two_port = TwoPort(; Phi = 0.0)
    @parameters begin
        rho = 0.098e-6
        l = 1
        A = 1
        R = rho * l / A # Electrical resistance of eddy current path
    end
    @equations begin
        D(Phi) ~ V_m * R
    end
end
EddyCurrent.f(; Phi, name, kwargs...) = EddyCurrent.f(; two_port__Phi = Phi, name, kwargs...)
