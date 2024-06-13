using Euler2D
using LinearAlgebra
using ShockwaveProperties
using Unitful
using ShockwaveDetection

ρL = 1.0u"kg/m^3"
vL = [0.0u"m/s"]
PL = 10.0u"Pa"
TL = uconvert(u"K", PL / (ρL * DRY_AIR.R))
ML = vL/speed_of_sound(ρL, PL; gas=DRY_AIR)

ρR = 0.125 * ρL
vR = [0.0u"m/s"]
PR = 0.1 * PL
TR = uconvert(u"K", PR/(ρR * DRY_AIR.R))
MR = vR/speed_of_sound(ρR, PR; gas=DRY_AIR)

s_high = PrimitiveProps(ρL, ML, TL)
s_low = PrimitiveProps(ρR, MR, TR)

sod1_bcs = EdgeBoundary(FixedPhantomOutside(s_high, DRY_AIR), FixedPhantomOutside(s_low, DRY_AIR))
copy_bcs = EdgeBoundary(ExtrapolateToPhantom(), ExtrapolateToPhantom())
u0_sod1(x) = ConservedProps(x < 0.5 ? s_high : s_low; gas=DRY_AIR) |> state_to_vector

flow_data = simulate_1D(0.0, 2.0, 2000, copy_bcs, 0.05, u0_sod1, DRY_AIR, 0.75)

println(flow_data)