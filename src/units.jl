# mass:        g / mol
# length:      nm
# time:        ps
# momentum:    g nm / ps mol
# energy:      kJ / mol
# temperature: K

const to_kilo = 1e-3
const to_pico = 1e12

# 2018 CODATA recommended values from NIST.
const N_A = 6.022_140_76e23 # 1 / mol
const planck = 6.626_070_15e-34 * to_kilo * to_pico * N_A # kJ ps / mol
const hbar = planck / 2pi
const kB = 1.380_649e-23 * to_kilo * N_A # kJ / mol K
