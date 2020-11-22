struct StepSize1D1 <: StepSize
end

Base.:*(ss::StepSize1D1, k::Float64) = ss


struct StepSize1Dj <: StepSize
    dq::Float64
end

Base.:*(ss::StepSize1Dj, k::Float64) = StepSize1Dj(ss.dq * k)


struct State1D{P} <: State{P}
    qs::Vector{Float64}
end

State1D(P::Int, xi_init::Float64) = State1D{P}(xi_init * ones(P))

Base.copy(st::State1D{P}) where {P} = State1D{P}(copy(st.qs))

path_xi(st::State1D) = st.qs


function density_diff(sys::PimcSystem, st1::State1D{P}, st2::State1D{P}, j::Int) where{P}
    pre_density = (sys.pot(st2.qs[j]) - sys.pot(st1.qs[j])) / P
    k = sys.mass * P / (2 * hbar^2 * sys.beta^2)

    for j_p in [bead_prev(P, j), bead_next(P, j)]
        link1 = (st1.qs[j] - st1.qs[j_p])^2
        link2 = (st2.qs[j] - st2.qs[j_p])^2
        pre_density += k * (link2 - link1)
    end

    pre_density
end

step_1(::PimcSystem, st::State1D, ss::StepSize1D1) = true, st

function step_j(sys::PimcSystem, st::State1D, j::Int, ss::StepSize1Dj)
    st_new = copy(st)
    st_new.qs[j] += ss.dq * randn()

    pre_density = density_diff(sys, st, st_new, j)

    if rand() < exp(-sys.beta * pre_density)
        true, st_new
    else
        false, st
    end
end


function deriv_estimator_pre(::State1D)
    dlogJ = 0.0
    chvar = [1.0]

    dlogJ, chvar
end

function deriv_estimator_1(sys::PimcSystem, chvar::Vector{Float64}, st::State1D, j::Int)
    dot(chvar, sys.force(st.qs[j]))
end

function deriv_estimator_2(chvar::Vector{Float64}, st::State1D{P}) where {P}
    dot(chvar, 2 .* st.qs[1] .- st.qs[2] .- st.qs[P])
end
