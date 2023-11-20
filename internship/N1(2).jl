#Noise 1 (Synaptic noise) in tonic firing, with high σ, low rate
using Revise
using DrWatson
using Logging

@quickactivate "SNNExamples"

using Distributions, Random 
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
Is = 10:5:25 #in pA
N_neurons = 5
E = SNN.AdEx(; N = N_neurons, param = SNN.AdExParameter(; El = El, τm = τm, Vt = Vt, Vr = Vr))
in_rate = 60Hz
σ = 1

D1 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate))
D2 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate))
50:500
rfn = (n)-> rand(Normal(1,0),n)
EE1 = SNN.SpikingSynapse(D1,E, :ge; rfn=rfn, p = 1.)
EE2 = SNN.SpikingSynapse(D2,E, :ge; rfn=rfn, p = 1.)
C = [EE1, EE2]
EE2.W

Ws = range(5,50,N_neurons) .* σ
@assert length(Ws) == length(EE1.W) 
for i in eachindex(Ws)
     EE1.W[i] = Ws[i]
     @info EE1.W[i]
end

heatmap(reshape(EE1.W, N_neurons, 1))
P = [D1, E]
SNN.monitor(E, [:v, :fire])

Random.seed!(11)
SNN.sim!(P, C; duration = 400ms)
v = SNN.getrecord(E, :v)
y = hcat(v...)'
mean(y, dims=1)
std(y,dims=1)
x = 1:length(v)
p = plot(x, y,
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "AdEx Neuron Model with Synaptic Noise \n (High σ, low rate)", 
leg = :none,
xaxis=("Time (ms)", extrema(x)),
yaxis=("Membrane Potential (mV)", extrema(y)),
palette= (:roma10) ,
lw = 4)