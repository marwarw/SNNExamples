using Revise
using DrWatson

@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots

default(
    size=(900, 600),
    tickfontsize=13,
    guidefontsize=18,
    margin=10Plots.mm,
    titlefontsize=18,
)
SNN.@load_units

Is = 1:10:200

const pA = ampere / 1e12

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV))

E.I = Is
E.v .= -60
SNN.monitor(E, [:v, :fire])

firing_rate = []
simtime = 1000ms
SNN.sim!([E], []; duration = simtime)

fr = sum(hcat(E.records[:fire]...),dims=2)

plot!(Is,fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "f-I Curve of AdEx Neuron",legend = :bottomright, label = "f-I Curve")


# Create scatter plot
scatter(Is, fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
    title = "Scatter Plot of Tonic Firing Response",
    legend = :bottomright, markersize = 4, color = :plum2, label = "Data Points")

