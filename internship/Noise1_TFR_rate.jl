#Noise 1 (Synaptic noise) in tonic firing, with high σ, high rate
using Revise
using DrWatson
using Logging

@quickactivate "SNNExamples"

using Distributions 
using SpikingNeuralNetworks
using Plots
using ColorSchemes
default(
    size=(900, 600),
    tickfontsize=12,
    guidefontsize=14,
    margin=10Plots.mm,
    titlefontsize=16,
)

SNN.@load_units

τm = 30ms
Vt = -50mV
Vr = -70mV
El = Vr
Is = 10:5:25 # in pA 
N_neurons = 5
E = SNN.AdEx(; N = N_neurons, param = SNN.AdExParameter(; El = El, τm = τm, Vt = Vt, Vr = Vr))
E.v .=El
in_rate = 200Hz 
σ = 0.
μ = 100.
rfn = (n)-> rand(Normal(μ,σ),n)
rfn(10)
P,C = [],[]
push!(P,E)
in_rates = range(5,50,N_neurons)*20
# @assert length(Ws) == length(EE1.W) 
for i in 1:N_neurons
    in_rate = in_rates[i]Hz
    D = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate))
    EE = SNN.SpikingSynapse(D,E, :ge; rfn=rfn,p = 1.)
    EE.W .=0; EE.W[i] = 10
    @info i, EE.W
    push!(C, EE)    
    push!(P, D)
end
length(P)
SNN.monitor(E, [:v, :fire])
SNN.sim!(P, C; duration = 1000ms)
v = SNN.getrecord(E, :v)
y = hcat(v...)'
mean(y, dims=1)
std(y,dims=1)
x = 1:length(v)
p = plot(x, y,
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "AdEx Neuron Model with Synaptic Noise \n (High σ, high rate)", 
leg = :none,
xaxis=("Time (ms)", extrema(x)),
yaxis=("Membrane Potential (mV)", extrema(y)),
palette= (:vanimo10) ,
lw = 4)

