abstract type CircEl end


struct Inductor <: CircEl
        ESR::Float64
        Induct::Float64
        ParCap::Float64
        Impedance::Float64
        Node1
        Node2
end

struct Capacitor <: CircEl
        ESR::Float64
        ESL::Float64
        Cap::Float64
        Impedance::Float64
        Node1
        Node2
end

struct Resistor <: CircEl
        R::Float64
        Impedance::Float64
        Node1
        Node2
        #NodePair(Node1,Node2) = new(Node1,Node2)
end
