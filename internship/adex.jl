#Spiking Neural Network simulation with AdEx neurons and Poisson inputs
using DrWatson
@quickactivate "SNNExamples"
using Plots
using SpikingNeuralNetworks

SNN.@load_units 

const pA = ampere / 1e12

I = SNN.Poisson(;N=1000,param=SNN.PoissonParameter(;rate=1Hz))
E = SNN.AdEx(;N=10,param=SNN.AdExParameter(; El=-70mV, Vt=-50.4mV, τm=20ms, a=2nS))
EE = SNN.SpikingSynapse(I,E, :ge; σ=50, p=0.2)
P=[I, E]

using Logging
 for x in fieldnames(SNN.AdExParameter)
     println(x, "  ",  getfield(E.param, x))
end

SNN.monitor(E, [:v, :fire, :w])
SNN.sim!(P, [EE]; duration = 700ms)
p1= plot(SNN.vecplot(E,:w),
xlabel = "Time (ms)", 
ylabel = "Adaptation current (pA)", 
title = "SNN Simulation with AdEx neurons and Poisson inputs")

p2= plot(SNN.vecplot(E,:v),
xlabel = "Time (ms)",
ylabel = "Membrane Potential (mV)", 
title = "SNN Simulation with AdEx neurons and Poisson inputs")


a = SNN.AdExParameter(;El=-49mV)
hline!([a.Vt],color=:red, label="Vt")
hline!([a.El],color=:blue, label="El")
plot!(ylims=(-70mV,-40mV))
E.v[1]=-70

plot((SNN.raster([E])); xlabel = "Time (ms)", ylabel = "AdEx neuron", title = "AdEx Neuron Spike Times", markersize=4, color=:oslo10, xlims = (0, 30000))


plot(SNN.vecplot(E,:v), 
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "SNN Simulation with AdEx neurons and Poisson inputs")