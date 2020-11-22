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
    "--l-max"
        metavar = "L"
        help = "basis truncation"
        arg_type = Int
        required = true
    "--dr"
        metavar = "D"
        help = "grid size"
        arg_type = Float64
        required = true
    "--grid-max"
        metavar = "G"
        help = "grid truncation"
        arg_type = Float64
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
    "--rho-file"
        metavar = "PATH"
        help = "path to density output file"
        arg_type = String
end
c = parse_args(ARGS, s, as_symbols=true)

const mass = c[:mass]
const epsilon = c[:epsilon]
const sigma = c[:sigma]
const l_max = c[:l_max]
const dr = c[:dr]
const grid_max = c[:grid_max]

if isnothing(c[:temp]) == isnothing(c[:temp_P])
    error("One of --temp or --temp-P must be specified.")
elseif !isnothing(c[:temp])
    const temperature = c[:temp]
else
    const temperature = c[:temp_P] / c[:P]
end

const P = c[:P]

const rho_file = c[:rho_file]

pot(r) = 4.0 * epsilon * ((sigma / r)^12 - (sigma / r)^6)

sys = NmmSystem(pot, mass, temperature)
np = NmmParameters(l_max, dr, grid_max)
grid, rho_diag = nmm(sys, np, P)

if !isnothing(rho_file)
    open(rho_file, "w") do f
        for i in eachindex(grid, rho_diag)
            println(f, "$(grid[i]) $(rho_diag[i])")
        end
    end
end
