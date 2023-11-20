using Revise
using DrWatson
@quickactivate "SNNExamples"
using Statistics
using SpikingNeuralNetworks
using Plots
using Interact
default(
    size=(900, 600),
    tickfontsize=10,
    guidefontsize=10,
    margin=8Plots.mm,
    titlefontsize=17,
    titlefontcolor=:orchid4,
)

SNN.@load_units

τm = 30ms   #membrane capacitance
Vt = -50mV  #initial threshold 
Vr = -70mV #reference resting potential
El = Vr  #resting membrane potential

Is = 10:5:25
Is  

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = τm, Vt = Vt, Vr = Vr))

SNN.monitor(E, [:v,:fire,:w])
E.v .= -70mV 
E.I = Is
my_Is = []
for x in 1:400
    E.I = E.I.+randn(4)
    SNN.sim!([E], []; duration = 1ms)
    push!(my_Is, E.I)
end
E.I .= 0

my_Is = hcat(my_Is...)

my_Is
length(1:10:4000)

q = plot(x->0,1:2000, label="No current")
[plot!(x->Is[i],2000:6000, label=string(Is[i]) , legend = true) for i in eachindex(Is)]

plot!(x->0,6000:10000, label="No current", xlabel="Time (ms)", ylabel="Input current (pA)", title="Plot input currents")

q2= plot(1:10:4000,my_Is[:,:]', legend= true, label= ["10pA" "15pA" "20pA" "25pA"], xlabel="Time (ms)", ylabel="Input currents (pA)", title="Input currents over time")

p= plot(SNN.vecplot(E,:v), legend= true, xlabel="Time (ms)", ylabel="Membrane potential (mV)", title="AdEx neuron in response to stochastic (noisy) input currents")

f= (plot(p, q2, layout=(2,1)))
plot!(; legend=:topright)

a_c= plot(SNN.vecplot(E,:w),
palette= (:vanimo10),
markersize=4,
color = :plum2,
xlabel = "Time (ms)", 
ylabel = "Adaptation current (pA)", 
title = "SNN Simulation with AdEx neurons ")
plot!(; legend=:topright)



w = hcat(E.records[:w]...)

avg_w = []
for i in 1:length(Is)
    w_current = E.records[:w][i, :]
    avg_w_current = mean(w_current)
    push!(avg_w, avg_w_current)
end

plot(Is, avg_w, markersize=10, markercolor=:blue, legend=true, xlabel="Input Current (pA)", ylabel="Average Adaptation Current", title="Average Adaptation Current vs. Input Current")


using Revise
using DrWatson
@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots
import Interact
using Interact
default(
    size=(600, 400),
    tickfontsize=8,
    guidefontsize=10,
    margin=5Plots.mm,
    titlefontsize=13,
    titlefontcolor=:teal,
    markersize=6,
    markercolor=:red,
)

SNN.@load_units 

# Default parameters
τm = 30ms  # Membrane capacitance
Vt = -50mV  # Initial threshold
Vr = -70mV  # Reference resting potential
El = Vr  # Resting membrane potential
Is = 10:5:25  # Input currents
initial_w = -65mV  # Initial adaptation current

# Create AdEx neuron model
E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = El, τm = τm, Vt = Vt, Vr = Vr))

# Monitor membrane potential and firing events
SNN.monitor(E, [:v, :fire, :w])

# Initialize membrane potential
E.v .= -70mV

# Set input currents
E.I = Is

# Simulate neuron response
SNN.sim!([E], []; duration = 1000ms)

# Extract membrane potential and adaptation current data
v = hcat(E.records[:v]...)
w = hcat(E.records[:w]...)

# Create interactive plot
plot((1:1000:4000, v, xlabel = "Time (ms)", ylabel = "Membrane potential (mV)", title = "Membrane potential vs. time"), w = -65:1:40, label = "Adaptation current")


