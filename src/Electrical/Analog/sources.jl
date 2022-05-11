"""
    Voltage(;name)   

Source for voltage,

"""
function Voltage(;name)   
    @named oneport = OnePort()
    @unpack v, i = oneport
    @named V = RealInput()
    eqs = [
        v ~ V.u
    ]
    
    extend(ODESystem(eqs, t, [], []; name=name, systems=[V]), oneport)
end
# Current Sources ######################################################################################################
"""
    Current(;name)   

Source for current.

"""
function Current(;name)   
    @named oneport = OnePort()
    @unpack v, i = oneport
    @named I = RealInput()
    eqs = [
        i ~ I.u
    ]
    
    extend(ODESystem(eqs, t, [], []; name=name, systems=[I]), oneport)
end
