#!/usr/bin/env julia

using SnaggedNecklace

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table! s begin
    "--mass"
        metavar = "M"
        help = "mass (g/mol)"
        arg_type = Float64
        required = true
    "--epsilon"
        metavar = "E"
        help = "LJ epsilon (kJ/mol)"
        arg_type = Float64
        required = true
    "--sigma"
        metavar = "S"
        help = "LJ sigma (nm)"
        arg_type = Float64
        required = true
    "--steps-equil"
        metavar = "N"
        help = "number of equilibration steps (default: 10% of production)"
        arg_type = Int
    "--steps-production"
        metavar = "N"
        help = "number of production steps"
        arg_type = Int
        required = true
    "--temp"
        metavar = "T"
        help = "temperature"
        arg_type = Float64
    "--temp-P"
        metavar = "TP"
        help = "temperature multiplied by P"
        arg_type = Float64
    "-P"
        metavar = "P"
        help = "number of beads"
        arg_type = Int
        required = true
    "--xi"
        metavar = "X"
        help = "value of constrained coordinate"
        arg_type = Float64
        required = true
    "--output-file"
        metavar = "PATH"
        help = "path to main output file"
        arg_type = String
    "--hist-file-template"
        metavar = "PATH"
        help = "template for paths to histogram output files"
        arg_type = String
end
c = parse_args(ARGS, s, as_symbols=true)

const mass = c[:mass]
const epsilon = c[:epsilon]
const sigma = c[:sigma]

if !isnothing(c[:steps_equil])
    const steps_equil = c[:steps_equil]
else
    const steps_equil = div(c[:steps_production], 10)
end

const steps_production = c[:steps_production]

if isnothing(c[:temp]) == isnothing(c[:temp_P])
    error("One of --temp or --temp-P must be specified.")
elseif !isnothing(c[:temp])
    const temperature = c[:temp]
else
    const temperature = c[:temp_P] / c[:P]
end

const P = c[:P]
const xi_init = c[:xi]

const output_file = c[:output_file]
const hist_file_template = c[:hist_file_template]

function pot(xs)
    r = sqrt(sum(xs.^2))
    4.0 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
end

function force(xyz::Vector{Float64})
    r2 = sum(xyz.^2)
    r = sqrt(r2)
    24.0 * epsilon * (2 * (sigma / r)^12 - (sigma / r)^6) .* xyz ./ r2
end

const step_size_1 = StepSize3D1(0.01, 0.01pi)
const step_size_j = StepSize3Dj(0.01)

sys = PimcSystem(pot, force, mass, temperature)
ip = PimcParameters(steps_equil, steps_production, step_size_1, step_size_j)
st = State3D(P, xi_init)

((mean1, err1, acsteps1), (mean2, err2, acsteps2), accept_ratio_1,
    accept_ratio_j, hists) = pimc(sys, ip, st)

if !isnothing(output_file)
    open(output_file, "w") do f
        print(f, "$(mean1) $(err1) $(acsteps1) ")
        print(f, "$(mean2) $(err2) $(acsteps2) ")
        println(f, "$(accept_ratio_1) $(accept_ratio_j)")
    end
end

if !isnothing(hist_file_template)
    for (j, (xis, ws)) in enumerate(hists)
        hist_file = replace(hist_file_template, "{}" => j)
        open(hist_file, "w") do f
            for (xi, w) in zip(xis, ws)
                println(f, "$(xi) $(w)")
            end
        end
    end
end
