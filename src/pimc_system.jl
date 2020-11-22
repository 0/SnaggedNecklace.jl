abstract type StepSize end

abstract type State{P} end

struct PimcSystem
    pot::Function
    force::Function
    mass::Float64
    beta::Float64

    function PimcSystem(pot::Function, force::Function, mass::Float64, temperature::Float64)
        beta = 1 / (kB * temperature)
        new(pot, force, mass, beta)
    end
end


bead_prev(P::Int, j::Int) = mod(j - 1, 1:P)
bead_next(P::Int, j::Int) = mod(j + 1, 1:P)
