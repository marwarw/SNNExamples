#Tonic firing neural response for increasing current step values
using Revise
using DrWatson

@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using ColorSchemes
using Plots
default(
    size=(900, 600),
    tickfontsize=13,
    guidefontsize=18,
    margin=10Plots.mm,
    titlefontsize=16,
)

SNN.@load_units 

#Define parameters
Is = [10, 20, 30, 40, 50]  # increasing current step values
firing_rate = []  #firing rate array

#Define AdEx neuron
E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El=-70mV, τm=20ms, τe=30ms, a=0.0nS, b=60pA, Vr=-55mV, Vt=-50mV, ΔT=2mV))
E.I = Is
SNN.monitor(E, [:v, :fire])

#Simulate and collect firing rates
for i in 1:length(Is)
    D = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = Is[i]))
    EE = SNN.SpikingSynapse(D, E, :ge; σ = 50, p = 1.0)
    P = [D, E]
    C = [EE]
    
    SNN.sim!(P, C; duration = 500ms)
    
    total_spikes = sum(E.fire)
    firing_rate_avg = total_spikes / 0.5
    push!(firing_rate, firing_rate_avg)
end

#Plot membrane potential over time
plot(SNN.vecplot(E, :v), 
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "Tonic firing neural response for increasing current step values",
leg = :none,
palette= (:berlin10))


# Create scatter plot
scatter(Is, firing_rate, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "Scatter Plot of Tonic Firing Response",
legend = false, markersize = 8, color = :blue)

##code

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

Is = 1:5:200

const pA = ampere / 1e12

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV))
E.I = Is
E.v .= -60
SNN.monitor(E, [:v, :fire])
firing_rate = []
simtime = 1000ms
SNN.sim!([E], []; duration = simtime)

fr = sum(hcat(E.records[:fire]...),dims=2)

p=plot!(Is,fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "f-I Curve of AdEx Neuron",legend = :bottomright, label = "f-I Curve")

s=scatter(Is, fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "Scatter Plot of Tonic Firing Response",
legend = :bottomright, markersize = 4, color = :plum2, label = "Data Points")

plot(p,s)

# Set up the AdEx neuron with initial parameters
const pA = ampere / 1e12
E = SNN.AdEx(; N = 1, param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV, a=0.01pA))

# Add the adaptation current manually
E.w .= 0.01pA

# Set up the range of adaptation current values
adapt_values = 0.01:0.01:0.1  # Adjust the range based on your requirements
E = SNN.AdEx(; N = 1, param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV, a=0.01pA))

SNN.monitor(E, [:v, :fire, :w])
E.w .=[]
SNN.sim!(E, []; duration = 1000ms)

p1= plot(SNN.vecplot(E,:w),
xlabel = "Time (ms)", 
ylabel = "Adaptation current (pA)", 
title = "Adex w")


# Initialize a plot
plot_handle = plot()

# Loop through each adaptation current value
for a_value in adaptation_current_values
    # Set the adaptation current value
    E.param.a .=[]

    # Simulate the AdEx neuron
    SNN.sim!([E], []; duration = 1000ms)

    # Extract the membrane potential data
    membrane_potential_data = hcat(E.records[:v]...)

    # Plot the membrane potential over time
    plot!(plot_handle, SNN.get_time_vector(E), membrane_potential_data, label="a = $a_value pA")
end

# Show the final plot
plot!(plot_handle, xlabel="Time (ms)", ylabel="Membrane Potential (mV)", legend=:topright)
display(plot_handle)



#TF_Is>

using Revise
using DrWatson

@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots

default(
    size=(900, 600),
    tickfontsize=10,
    guidefontsize=12,
    margin=10Plots.mm,
    titlefontsize=18,
    titlefontcolor=:teal,
)
SNN.@load_units

Is = 1:5:200

const pA = ampere / 1e12

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV))
E.I = Is
E.v .= -60
SNN.monitor(E, [:v, :fire, :w])

firing_rate = []
simtime = 1000ms
SNN.sim!([E], []; duration = simtime)

fr = sum(hcat(E.records[:fire]...),dims=2)

p=plot(Is,fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "f-I Curve of AdEx Neuron",legend = :bottomright, label = "f-I Curve")

s=scatter(Is, fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "Scatter Plot of Tonic Firing Response",
legend = :bottomright, markersize = 3, color = :plum2, label = "Data Points")

q= plot(p,s, layout=(2, 1))

v= hcat(E.records[:v]...)
V_subsampled = E.records[:v][1:1000:end]
mp= plot(1:100:4000, V_subsampled, legend= true, xlabel="Time (ms)",
ylabel="mV", title="MP over time")



#+adaptation_current

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

Is = 1:5:500

const pA = ampere / 1e12

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV))
E.I = Is
E.v .= -60
SNN.monitor(E, [:v, :fire])

firing_rate = []
simtime = 1000ms
SNN.sim!([E], []; duration = simtime)

fr = sum(hcat(E.records[:fire]...), dims=2)

fr_values = [];
w_values = [];

for a in 0:0.1:2
    param = SNN.AdExParameter(; El=-70mV, τm=20ms, Vt=-50mV, Vr=-70mV, a=a)
    E.param = param
    SNN.sim!([E], []; duration = simtime)
    fr = sum(hcat(E.records[:fire]...))
    println(a, " ", fr)
    push!(fr_values, fr)
    push!(w_values,a)
end



p= plot(w_values, fr_values, xlabel = "Adaptation Current (w)", ylabel = "Firing Rate (Hz)", title = "f-w Curve of AdEx Neuron")

for b in 0:2:60
    param = SNN.AdExParameter(; El=-70mV, τm=20ms, Vt=-50mV, Vr=-70mV, a=2nS, b=b)
    E.param = param
    SNN.sim!([E], []; duration = simtime)
    fr = sum(hcat(E.records[:fire]...))
    push!(fr_values, fr)
    push!(w_values,b)
end

q= plot!(w_values, fr_values, xlabel = "Adaptation Current (w)", ylabel = "Firing Rate (Hz)", title = "f-w Curve of AdEx Neuron")

m= plot(SNN.vecplot(E,:v), 
xlabel = "Time (ms)", 
ylabel = "Membrane Potential (mV)", 
title = "Membrane Potential vs. Time")

#

using Revise
using DrWatson

@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots

default(
  size=(900, 600),
  tickfontsize=10,
  guidefontsize=12,
  margin=10Plots.mm,
  titlefontsize=18,
  titlefontcolor=:teal,
)
SNN.@load_units

Is = 1:5:200

const pA = ampere / 1e12

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV))
E.I = Is
E.v .= -60

SNN.monitor(E, [:v, :fire])

firing_rate = []
simtime = 2000ms
dt = 50ms  # Sampling time step
steps = round(simtime / dt)  # Number of time steps
time= 1:dt:simtime  # Time vector

SNN.sim!([E], []; duration = simtime)

fr = sum(hcat(E.records[:fire]...),dims=2)

V= sum(hcat(E.records[:v]...), dims=2) # Extract membrane potential data

p=plot(Is,fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "f-I Curve of AdEx Neuron",legend = :bottomright, label = "f-I Curve")

s=scatter(Is, fr, xlabel = "Input Current (pA)", ylabel = "Firing Rate (Hz)",
title = "Scatter Plot of Tonic Firing Response",
legend = :bottomright, markersize = 3, color = :plum2, label = "Data Points")

q= plot(p,s, layout=(2, 1))


m= plot(time, V, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)",
title = "Membrane Potential vs. Time", legend = :bottomright, label = "Membrane Potential")
