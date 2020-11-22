struct NmmSystem
    pot::Function
    mass::Float64
    beta::Float64

    function NmmSystem(pot::Function, mass::Float64, temperature::Float64)
        beta = 1 / (kB * temperature)
        new(pot, mass, beta)
    end
end
