struct NmmParameters
    l_max::Int
    dr::Float64
    grid_max::Float64
end

function nmm(sys::NmmSystem, np::NmmParameters, P::Int)
    grid = collect(range(np.dr, np.grid_max; step=np.dr))

    eK_m = zeros(length(grid), length(grid))
    eK_p = zeros(length(grid), length(grid))
    k = sys.mass * P / (2 * hbar^2 * sys.beta)
    for i in eachindex(grid)
        for j in eachindex(grid)
            eK_m[j, i] += exp(-k * (grid[i] - grid[j])^2)
            eK_p[j, i] += exp(-k * (grid[i] + grid[j])^2)
        end
    end
    eK = eK_m - eK_p
    eK ./= maximum(abs.(eK))

    scaling = 0.0
    scaling_event = Threads.Event()

    rho_diag = zeros(length(grid))
    add_lock = Threads.SpinLock()

    Threads.@threads for l in 0:np.l_max
        eV = zeros(length(grid), length(grid))
        for i in eachindex(grid)
            pot_rot = hbar^2 * l * (l + 1) / (2 * sys.mass * grid[i]^2)
            eV[i, i] += exp(-0.5 * sys.beta / P * (sys.pot(grid[i]) + pot_rot))
        end

        eH = eV * eK * eV
        F = eigen(Symmetric(eH))

        if l == 0
            scaling = maximum(abs.(F.values))
            notify(scaling_event)
        else
            wait(scaling_event)
        end

        rho_diag_l = zeros(length(F.values))
        for i in eachindex(rho_diag_l)
            for j in eachindex(rho_diag_l)
                rho_diag_l[i] += F.vectors[i, j]^2 * (F.values[j] / scaling)^P
            end
        end

        lock(add_lock) do
            rho_diag .+= (2l+1) .* rho_diag_l
        end
    end

    rho_diag ./= sum(rho_diag) * (grid[2] - grid[1])

    grid, rho_diag
end
