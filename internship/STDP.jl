#Connectivity reflects coding: a model of voltage-based STDP with homeostasis
#Claudia Clopath et al., 2009
"https://www.nature.com/articles/nn.2479"

using Revise
using DrWatson
using Logging

@quickactivate "SNNExamples"

using Distributions 
using SpikingNeuralNetworks
using Plots
using ColorSchemes

SNN.@load_units

#Parameters
C= 281pF #membrane capacitance
gL= 30nS #leak conductance
El= -70.6mV #resting potential
ΔT= 2mV #slope factor
Vt= -50.4mV #threshold potential at rest
τw= 144ms #adaptation time constant
a= 4nS #subthreshold adaptation
b= 0.805pA #spike triggered adaptation
τm= 50ms #tVT, threshold potential time constant
#spike current after a spike
Vr= 30.4mV #VTmax, threshold potential after a spike
τm= 20ms #C/gL
τe= 5ms
τi= 10ms
Vr= -70.6mV

R= 0.500 #GΩ

function VoTri(rho, Dt, par)
    # Parameters of the plasticity rule
    A_m = par[1]                    # Amplitude for depression
    A_p = par[2]                    # Amplitude for potentiation
    tau_p = par[3]                  # Time constant for voltage filter for potentiation
    tetam = -70.6                   # Low voltage threshold
    tetap = -45.3                   # High voltage threshold
    tau_r = par[4]                  # Time constant low pass r [ms]
    tau_d = par[5]                  # Low pass of the membrane pot [ms] for the depression term

    # Protocol of Sjoestroem et al. Neuron 2001: pre and post spike trains
    n = 5                           # Number of pairings
    f = 1000 / rho                  # Time between two pairs [ms]
    l = n * f + abs(Dt) + 1         # Length of the spike trains
    x = zeros(Float64, l)           # Presynaptic spike train
    x[abs(Dt) + 1: f: end-1] .= 1
    y = zeros(Float64, l)           # Postsynaptic spike train
    y[abs(Dt) + 1 + Dt: f: end] .= 1

    # Protocol: repetition at 0.1Hz
    rho_low = 0.1
    n_rep = 15
    if rho == 0.1
        n_rep = 10
    end
    x = repeat([x, zeros(Float64, 1000 ÷ rho_low - length(x))], n_rep)
    y = repeat([y, zeros(Float64, 1000 ÷ rho_low - length(y))], n_rep)

    # Parameters of the neuron model
    E_L = -70.6                     # [mV] Resting potential

    # Initialization
    l = length(x)                   # Simulation length
    I_ext = zeros(Float64, l)       # Extra current injected
    I_s = y .* 1000000              # Current injected for postsynaptic spike induction
    w_tail = zeros(Float64, l)      # Current for spike after-depolarization
    V_T = -50.4                     # Adaptive threshold
    I_tot = zeros(Float64, l)       # Total current
    wad = zeros(Float64, l)         # Adaptation
    u = fill(E_L, l)                # Membrane potential
    u_md = fill(E_L, l)             # Filtered membrane potential for the depression term
    u_mp = fill(E_L, l)             # Filtered membrane potential for the potentiation term
    r = zeros(Float64, size(x))     # Low-pass the presynaptic spike train, x
    w = 0.5                         # Initialization of weights
    counter = 0                     # Trick to count how long is the spike - here spikes are forced to be 2ms long

    # Main loop over time
    for t = 4:l
        I_tot[t] = I_s[t] + I_ext[t]  # Current
        u[t], wad[t], w_tail[t], counter, V_T = aEIF(u[t-1], wad[t-1], w_tail[t-1], I_tot[t], counter, V_T)  # Voltage with the aEIF neuron model
        u_md[t+1] = u[t] / tau_d + (1 - (1 / tau_d)) * u_md[t]  # Low-pass of voltage for depression
        u_mp[t+1] = u[t] / tau_p + (1 - (1 / tau_p)) * u_mp[t]  # Low-pass of voltage for potentiation
        r[t+1] = x[t] / tau_r + (1 - (1 / tau_r)) * r[t]  # Low-pass of presynaptic spike
        u_sig = (u[t] > tetap) * (u[t] - tetap)  # Voltage threshold
        u_md_sig = ((u_md[t-3] - tetam) > 0) * (u_md[t-3] - tetam)  # Low-pass threshold
        u_mp_sig = ((u_mp[t-3] - tetam) > 0) * (u_mp[t-3] - tetam)  # Low-pass threshold
        w = w - A_m * x[t] * u_md_sig + A_p * u_sig * r[t] * u_mp_sig  # Weight updates
    end

    return w
end

