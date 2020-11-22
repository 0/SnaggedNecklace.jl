struct StepSize3D1 <: StepSize
    dcostheta::Float64
    dphi::Float64
end

function Base.:*(ss::StepSize3D1, k::Float64)
    dcostheta_k = min(0.5, ss.dcostheta * k)
    dphi_k = min(0.5pi, ss.dphi * k)
    StepSize3D1(dcostheta_k, dphi_k)
end


struct StepSize3Dj <: StepSize
    dr::Float64
end

Base.:*(ss::StepSize3Dj, k::Float64) = StepSize3Dj(ss.dr * k)


struct State3D{P} <: State{P}
    xyzs::Matrix{Float64}
end

function State3D(P::Int, xi_init::Float64)
    xyzs = zeros(Float64, 3, P)
    xyzs[3, :] .= xi_init
    State3D{P}(xyzs)
end

Base.copy(st::State3D{P}) where {P} = State3D{P}(copy(st.xyzs))

path_xi(st::State3D{P}) where {P} = [cart2sphr(st.xyzs[:, j])[1] for j in 1:P]


function cart2sphr(xyz::Vector{Float64})
    r = sum(xyz.^2)^(1/2)
    costheta = xyz[3] / r
    phi = atan(xyz[2], xyz[1])
    [r, costheta, phi]
end

function sphr2cart(r::Float64, costheta::Float64, phi::Float64)
    x = r * sqrt(1 - costheta^2) * cos(phi)
    y = r * sqrt(1 - costheta^2) * sin(phi)
    z = r * costheta
    [x, y, z]
end


function density_diff(sys::PimcSystem, st1::State3D{P}, st2::State3D{P}, j::Int) where {P}
    pre_density = (sys.pot(st2.xyzs[:, j]) - sys.pot(st1.xyzs[:, j])) / P
    k = sys.mass * P / (2 * hbar^2 * sys.beta^2)

    for j_p in [bead_prev(P, j), bead_next(P, j)]
        link1 = (st1.xyzs[:, j] .- st1.xyzs[:, j_p]).^2 |> sum
        link2 = (st2.xyzs[:, j] .- st2.xyzs[:, j_p]).^2 |> sum
        pre_density += k * (link2 - link1)
    end

    pre_density
end

function step_1(sys::PimcSystem, st::State3D, ss::StepSize3D1)
    r, costheta, phi = cart2sphr(st.xyzs[:, 1])
    costheta += ss.dcostheta * (rand() - 0.5)
    phi += ss.dphi * randn()

    if costheta < -1
        costheta = -2 - costheta
    elseif costheta > 1
        costheta = 2 - costheta
    end

    st_new = copy(st)
    st_new.xyzs[:, 1] .= sphr2cart(r, costheta, phi)

    pre_density = density_diff(sys, st, st_new, 1)

    if rand() < exp(-sys.beta * pre_density)
        true, st_new
    else
        false, st
    end
end

function step_j(sys::PimcSystem, st::State3D, j::Int, ss::StepSize3Dj)
    step_costheta = 2 * (rand() - 0.5)
    step_phi = 2pi * rand()
    step_r = ss.dr * rand()
    st_new = copy(st)
    st_new.xyzs[:, j] .+= sphr2cart(step_r, step_costheta, step_phi)

    pre_density = density_diff(sys, st, st_new, j)

    if rand() < exp(-sys.beta * pre_density)
        true, st_new
    else
        false, st
    end
end


function deriv_estimator_pre(st::State3D)
    xyz1 = st.xyzs[:, 1]
    r1, _, _ = cart2sphr(xyz1)

    dlogJ = 2 / r1
    chvar = xyz1 ./ r1

    dlogJ, chvar
end

function deriv_estimator_1(sys::PimcSystem, chvar::Vector{Float64}, st::State3D, j::Int)
    dot(chvar, sys.force(st.xyzs[:, j]))
end

function deriv_estimator_2(chvar::Vector{Float64}, st::State3D{P}) where {P}
    dot(chvar, 2 .* st.xyzs[:, 1] .- st.xyzs[:, 2] .- st.xyzs[:, P])
end
