#Simulation of an AdEx neuron model with synaptic noise
#Noise 1 (Synaptic noise) in tonic firing, with low σ and high firing rate
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

#Tonic firing neuron model parameters
τm = 30ms
Vt = -50mV
Vr = -70mV
El = Vr

#Define current inputs range and n. neurons
Is = 10:5:25 #in pA
N_neurons = 5

#Define AdEx neuron
E = SNN.AdEx(; N = N_neurons, param = SNN.AdExParameter(; El = El, τm = τm, Vt = Vt, Vr = Vr))

##Low σ, high rate
#Low σ = narrow distribution of the synaptic noise 
# = less variability in the synaptic inputs received by neurons

#Create Poisson populations and connections
in_rate = 200Hz
σ = 20

D1 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate))
D2 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate))

50:500
rfn = (n)-> rand(Normal(10,0),n)
EE1 = SNN.SpikingSynapse(D1,E, :ge; rfn=rfn,  σ = σ, p = 1.)
EE2 = SNN.SpikingSynapse(D2,E, :ge; rfn=rfn,  σ = σ, p = 1.)

C = [EE1, EE2]

#Ensure length of weight values (Ws) = length of synaptic weights(EE1.W)
Ws = range(5,50,N_neurons)
@assert length(Ws) == length(EE1.W) 
for i in eachindex(Ws)
     EE1.W[i] = Ws[i]
     @info EE1.W[i]
end

#Visualize EE1.W
heatmap(reshape(EE1.W, N_neurons, 1))

#Combine the populations and simulate
P = [D1, E]
SNN.monitor(E, [:v, :fire])
SNN.sim!(P, C; duration = 400ms)

#Retrieve + process recorded membrane potentials
v = SNN.getrecord(E, :v)
y = hcat(v...)'
mean(y, dims=1)
std(y,dims=1)
x = 1:length(v)


#Plot membrane potential over time
p = plot(x, y,
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "AdEx Neuron Model with Synaptic Noise \n (Low σ, high rate)", 
leg = :none,
xaxis=("Time (ms)", extrema(x)),
yaxis=("Membrane Potential (mV)", extrema(y)),
palette= (:berlin10) ,
lw = 4)