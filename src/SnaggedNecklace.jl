module SnaggedNecklace

using LinearAlgebra: Symmetric, dot, eigen, normalize
using Statistics: mean, std
using StatsBase: Histogram, fit

export
    NmmSystem,
    NmmParameters,
    nmm,

    PimcSystem,
    StepSize1D1,
    StepSize1Dj,
    State1D,
    StepSize3D1,
    StepSize3Dj,
    State3D,
    PimcParameters,
    pimc

include("units.jl")
include("stats.jl")

include("nmm_system.jl")
include("nmm.jl")

include("pimc_system.jl")
include("pimc_1d.jl")
include("pimc_3d.jl")
include("pimc_integrate.jl")

end
