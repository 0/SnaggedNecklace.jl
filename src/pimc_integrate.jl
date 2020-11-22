struct PimcParameters
    steps_equil::Int
    steps_production::Int
    step_size_1::StepSize
    step_size_j::StepSize
end

function adjustment_factor(acceptance_ratio)
    if acceptance_ratio < 0.2
        0.9
    elseif 0.2 <= acceptance_ratio < 0.3
        0.95
    elseif 0.3 <= acceptance_ratio < 0.4
        0.99
    elseif 0.4 == acceptance_ratio
        1.0
    elseif 0.4 < acceptance_ratio <= 0.5
        1.01
    elseif 0.5 < acceptance_ratio <= 0.6
        1.05
    elseif 0.6 < acceptance_ratio
        1.1
    end
end

function pimc(sys::PimcSystem, ip::PimcParameters, st::State{P}) where {P}
    sss = [ip.step_size_1; [ip.step_size_j for _ in 2:P]]

    nums_total = [0 for _ in 1:P]
    nums_accept = [0 for _ in 1:P]

    for i in 1:ip.steps_equil
        if i % 100 == 0
            for j in 1:P
                sss[j] *= adjustment_factor(nums_accept[j] / nums_total[j])
                nums_total[j] = 0
                nums_accept[j] = 0
            end
        end

        accept, st = step_1(sys, st, sss[1])
        nums_total[1] += 1
        accept && (nums_accept[1] += 1)

        for j in 2:P
            accept, st = step_j(sys, st, j, sss[j])
            nums_total[j] += 1
            accept && (nums_accept[j] += 1)
        end
    end

    num_total_1 = 0
    num_total_j = 0
    num_accept_1 = 0
    num_accept_j = 0

    derivs1 = zeros(ip.steps_production)
    derivs2 = zeros(ip.steps_production)
    dist = zeros(P, ip.steps_production)

    for i in 1:ip.steps_production
        accept, st = step_1(sys, st, sss[1])
        num_total_1 += 1
        accept && (num_accept_1 += 1)

        for j in 2:P
            accept, st = step_j(sys, st, j, sss[j])
            num_total_j += 1
            accept && (num_accept_j += 1)
        end

        dlogJ, chvar = deriv_estimator_pre(st)

        # Free energy derivative estimator 1.
        derivs1[i] += dlogJ

        for j in 1:P
            derivs1[i] += deriv_estimator_1(sys, chvar, st, j) * sys.beta / P
        end

        # Free energy derivative estimator 2.
        derivs2[i] += dlogJ
        derivs2[i] += deriv_estimator_1(sys, chvar, st, 1) * sys.beta / P
        derivs2[i] -= deriv_estimator_2(chvar, st) * sys.mass * P / (hbar^2 * sys.beta)

        dist[:, i] .= path_xi(st)
    end

    hists = NTuple{2,Vector{Float64}}[]

    for j in 1:P
        h = normalize(fit(Histogram, dist[j, :], nbins=50))
        edges = only(h.edges)
        mids = (edges[begin:(end-1)] .+ edges[(begin+1):end]) ./ 2
        push!(hists, (mids, h.weights))
    end

    (aggregate(derivs1), aggregate(derivs2), num_accept_1 / num_total_1,
        num_accept_j / num_total_j, hists)
end
