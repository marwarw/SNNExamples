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
using Statistics

const pA = ampere / 1e12
duration= 2000 #2s
Is = 10:10:100  # Increasing current step values

# AdEx neuron
function initialize_pop(firing_pattern)
    if firing_pattern == "tonic"
        param = SNN.AdExParameter(;El=-70mV, τm=20ms, τe=30ms, a=0.0nS, b=60pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "adapting"
        param = SNN.AdExParameter(;El=-70mV, τm=200ms, τe=100ms, a=0.0nS, b=5pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "bursting"
        param = SNN.AdExParameter(;El=-70mV, τm=5.0ms, τe=100ms, a=-0.5nS, b=7pA, Vr=-46mV, Vt=-50mV, ΔT=2mV)
    else
        error("Invalid firing pattern")
    end

    E = SNN.AdEx(;N=length(Is), param=param)
    SNN.monitor(E, [:v, :fire, :w])
    return E
end


Fp = ["tonic", "adapting", "bursting"]  # Different firing patterns to simulate
firing_rates = Dict{String, Vector{Float64}}()

# Tonic firing

E1 = initialize_pop("tonic")
firing_rates["tonic"] = Float64[]
for I in Is
    E1.I .= I
    SNN.sim!([E1], []; duration= 2000ms)
    rate = mean(E1.fire) / duration  # Calculate firing rate per 100ms
    push!(firing_rates["tonic"], rate)
end
p1 = plot(Is, firing_rates["tonic"], xlabel="Current (pA)", ylabel="Firing rate (Hz)", linecolor = :blue, label="Tonic firing")

# Adapting firing
E2 = initialize_pop("adapting")
firing_rates["adapting"] = Float64[]
for I in Is
    E2.I .= I
    SNN.sim!([E2], []; duration = 2000ms)
    rate = mean(E2.fire) / duration  # Calculate firing rate per 100ms
    push!(firing_rates["adapting"], rate)
end
p2 = plot(Is, firing_rates["adapting"], xlabel="Current (pA)", ylabel="Firing rate (Hz)", linecolor = :green, label="Adapting firing")

# Bursting firing
E3 = initialize_pop("bursting")
firing_rates["bursting"] = Float64[]
for I in Is
    E3.I .= I
    SNN.sim!([E3], []; duration = 2000ms)
    rate = mean(E3.fire) / duration  # Calculate firing rate per 100ms
    push!(firing_rates["bursting"], rate)
end
p3 = plot(Is, firing_rates["bursting"], xlabel="Current (pA)", ylabel="Firing rate (Hz)", linecolor = :red, label="Bursting firing")

E1.v .=E2.v .=E3.v .= -70mV

# Plot
plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))
plot!(title="Firing Rate vs Current Step for Different Firing Patterns")









##
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
using Statistics

const pA = ampere / 1e12

Is= 10:10:80  # Increasing current step values
Is

# AdEx neuron
function initialize_pop(firing_pattern)
    if firing_pattern == "tonic"
        param = SNN.AdExParameter(;El=-70mV, τm=20ms, τe=30ms, a=0.0nS, b=60pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "adapting"
        param = SNN.AdExParameter(;El=-70mV, τm=200ms, τe=100ms, a=0.0nS, b=5pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "bursting"
        param = SNN.AdExParameter(;El=-70mV, τm=5.0ms, τe=100ms, a=-0.5nS, b=7pA, Vr=-46mV, Vt=-50mV, ΔT=2mV)
    else
        error("Invalid firing pattern")
    end

    E = SNN.AdEx(;N=length(Is), param=param)
    SNN.monitor(E, [:v, :fire, :w])
    return E
end


Fp = ["tonic", "adapting", "bursting"]  # Different firing patterns to simulate
Firing_rates = Dict{String, Vector{Float64}}()

# Tonic firing
E1 = initialize_pop("tonic")
E1.I = Is
SNN.sim!([E1], []; duration = 100ms)
p1 =plot(SNN.vecplot(E1, :w), SNN.vecplot(E1, :v), label="Tonic firing")

# Adapting firing
E2 = initialize_pop("adapting")
E2.I = Is
SNN.sim!([E2], []; duration = 100ms)
p2 = plot(SNN.vecplot(E2, :w), SNN.vecplot(E2, :v), label="Adapting firing")

# Bursting firing
E3 = initialize_pop("bursting")
E3.I = Is
SNN.sim!([E3], []; duration = 100ms)
p3 = plot(SNN.vecplot(E3, :w), SNN.vecplot(E3, :v), label="Bursting firing")

#Plot
plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))
xlabel!("Current (pA)")
ylabel!("Firing rate (Hz)")
plot!(title="Multiple Firing Patterns \n in AdEx Neuron Model")






for pattern in firing_patterns
    firing_rates[pattern] = Float64[]
    for current_step in current_steps
        E = initialize_pop(pattern)
        
        I = SNN.Poisson(; N = length(current_steps), param = SNN.PoissonParameter(rate = current_step))
        synapse = SNN.SpikingSynapse(I, E, :ge; σ = 0.01, p = 1.0)
        P = [I, E]
        C = [synapse]
        
        SNN.sim!(P, C; duration = 1000ms)
        
        firing_rate_avg = mean(E.fire) / 1.0  # Divide by simulation duration in seconds
        push!(firing_rates[pattern], firing_rate_avg)
    end
end


# Tonic firing
E1 = initialize_pop("tonic")
E1.I = [40]
SNN.sim!([E1], []; duration = 500ms)
p1 = plot(SNN.vecplot(E1, :w), SNN.vecplot(E1, :v), label="Tonic firing")



p1= plot(current_steps, firing_rates["tonic"], marker = :circle, linecolor = :blue, label = "Tonic")
p2= plot!(current_steps, firing_rates["adapting"], marker = :circle, linecolor = :green, label = "Adapting")
p3= plot!(current_steps, firing_rates["bursting"], marker = :circle, linecolor = :purple, label = "Bursting")
plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))
xlabel!("Current Step (Is)")
ylabel!("Firing Rate (spikes/s)")
title!("Neural Response Firing Rate for Increasing Current Step")









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

# Function to initialize the population with specified firing pattern
function initialize_pop(firing_pattern)
    if firing_pattern == "tonic"
        param = SNN.AdExParameter(;El=-70mV, τm=20ms, τe=30ms, a=0.0nS, b=60pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "adapting"
        param = SNN.AdExParameter(;El=-70mV, τm=200ms, τe=100ms, a=0.0nS, b=5pA, Vr=-55mV, Vt=-50mV, ΔT=2mV)
    elseif firing_pattern == "bursting"
        param = SNN.AdExParameter(;El=-70mV, τm=5.0ms, τe=100ms, a=-0.5nS, b=7pA, Vr=-46mV, Vt=-50mV, ΔT=2mV)
    else
        error("Invalid firing pattern")
    end

    E = SNN.AdEx(;N=1, param=param)
    return E
end

# Simulate firing rate for increasing current step
function simulate_firing_rate(firing_pattern, current_steps, sim_duration)
    firing_rates = []
    
    for I in current_steps
        E = initialize_pop(firing_pattern)
        E.I = [I]
        
        I_pop = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(rate = I))
        synapse = SNN.SpikingSynapse(I_pop, E, :ge; σ = 0.01, p = 1.0)
        
        SNN.monitor(E, [:fire])
        SNN.sim!([I_pop, E], [synapse]; duration = 1000ms)
        
        total_spikes = sum(E.fire)
        firing_rate = total_spikes / (1000 / 1000)  # Divide by simulation duration in seconds
        push!(firing_rates, firing_rate)
    end
    
    return firing_rates
end

# Define the parameters
current_steps = [10, 20, 30, 40, 50]  # Increasing current step values

# Simulate firing rate for each firing pattern
tonic_firing_rates = simulate_firing_rate("tonic", current_steps, 500ms)
adapting_firing_rates = simulate_firing_rate("adapting", current_steps, 300ms)
bursting_firing_rates = simulate_firing_rate("bursting", current_steps, 30ms)

# Plot firing rate for each firing pattern
p1= plot!(current_steps, tonic_firing_rates, marker = :circle, xlabel = "Current Step (I)", ylabel = "Firing Rate (spikes/s)", label = "Tonic firing")
p2= plot!(current_steps, adapting_firing_rates, marker = :circle, label = "Adapting firing")
p3= plot!(current_steps, bursting_firing_rates, marker = :circle, label = "Bursting firing")
plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))


using Revise
using DrWatson
@quickactivate "SNNExamples"

using Plots
using SpikingNeuralNetworks
SNN.@load_units

# Define the parameters
Is = [10,30,50,70,90]  # Increasing current step values
firing_rate = []  #firing rate array
rate = [1,5,10,15,20]
rate_change = (20-1)/5

# Define the Poisson population with rate change
I= SNN.Poisson(;N = length(Is), param = SNN.PoissonParameter(;rate = rate_change))

# Define the AdEx neuron
E= SNN.AdEx(;N = length(Is), param = SNN.AdExParameter(;El = -70mV, τm = 30ms, Vt = -50mV, Vr = -70mV))
E.I= Is

# Create the SpikingSynapse connection
EE= SNN.SpikingSynapse(I, E, param = SNN.SpikingSynapseParameter(;Wmax = 0.01), :ge; σ = 0.01, p = 1.0)
#S.param = SNN.SpikingSynapseParameter(;Wmax = 0.01)

# Perform the simulation and collect firing rates
for i in 1:length(Is)
    I = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = Is[i]))
    EE = SNN.SpikingSynapse(I, E, :ge; σ = 0.01, p = 1.0)
    P = [I, E]
    C = [EE]
    
    SNN.sim!(P, C; duration = 1000ms)
    
    total_spikes = sum(E.fire)
    firing_rate_avg = total_spikes / 1.0  # Divide by simulation duration in seconds
    push!(rate, firing_rate_avg)
end

# Plot the firing rates
plot(Is, firing_rate, marker = :circle, xlabel = "Current Step (Is)", ylabel = "Firing Rate (spikes/s)", legend = false)

# Combine the populations and connections
P = [I, E]
C = [EE]

# Perform the simulation
SNN.monitor(E, [:v, :fire])
SNN.sim!(P, C; duration = 1000ms)
SNN.vecplot(E, :fire)


# Extract the firing rates of the AdEx population
firing_rates = SNN.getrecord(E, :fire)

# Plot the firing rate as a function of increasing step current
plot(Is, firing_rates, xlabel = "Increasing Step Current (Is)", ylabel = "Firing Rate")



# Plot the firing rate as a function of the input current
plot(input_currents, SNN.rate(neurons), marker = :circle, xlabel = "Input Current (A)", ylabel = "Firing Rate (Hz)", legend = false)





using Revise
using DrWatson
@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots
default(
    size=(900,600),
    tickfontsize=13,
    guidefontsize=18,
    margin=10Plots.mm,
    titlefontsize=18,
)
SNN.@load_units


# AdEx neuron
function initialize_pop(firing_pattern)
    if firing_pattern == "tonic"
        param = SNN.AdExParameter(;El=-70mV, τm = 20ms, τe = 30ms, a=0.0nS, b=0.06nA, Vr = -55mV)
    elseif firing_pattern == "adapting"
        param = SNN.AdExParameter(;El=-70mV, τm = 200ms, τe = 100ms, a=0.0nS, b=0.005nA, Vr = -55mV)
    elseif firing_pattern == "bursting"
        param = SNN.AdExParameter(;El=-70mV, τm = 5.0ms, τe = 100ms, a=0.5nS, b=0.007nA, Vr = -51mV)
    else
        error("Invalid firing pattern")
    end

    E = SNN.AdEx(; N=1, param=param)
    SNN.monitor(E, [:v, :fire, :w])
    return E
end

# Tonic firing
E1 = initialize_pop("tonic")
E1.I = [65]
SNN.sim!([E1], []; duration = 300ms)
p1 = plot(SNN.vecplot(E1, :w), SNN.vecplot(E1, :v), label="Tonic firing")

# Adapting firing
E2 = initialize_pop("adapting")
E2.I = [65]
SNN.sim!([E2], []; duration = 300ms)
p2 = plot(SNN.vecplot(E2, :w), SNN.vecplot(E2, :v), label="Adapting firing")

# Bursting firing
E3 = initialize_pop("bursting")
E3.I = [65,-65,65,-65,65,-65,65,-65]
SNN.sim!([E3], []; duration = 300ms)
p3 = plot(SNN.vecplot(E3, :w), SNN.vecplot(E3, :v), label="Bursting firing")

plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))
xlabel!("Time(ms)")
ylabel!("Membrane Potential (mV)")
plot!(title="Multiple Firing Patterns in AdEx Neuron Model")















#Multiple firing patterns in the AdEx neuron model
"https://neuronaldynamics.epfl.ch/online/Ch6.S1.html"

using Revise
using DrWatson
@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots
default(
    size=(900,600),
    tickfontsize=13,
    guidefontsize=18,
    margin=10Plots.mm,
    titlefontsize=18,
)
SNN.@load_units


# AdEx neuron
function initialize_pop()
    E = SNN.AdEx(;N = 1, param=SNN.AdExParameter(;El=-70mV, τm = 200ms, τe = 100ms, a=0.0nS, b=0.005nA, Vr = -55mV))
    SNN.monitor(E, [:v, :fire, :w])
    for param=SNN.AdExParameter(;El=-70)
        if SNN.AdExParameter(;El=-70mV, τm = 20ms, τe = 30ms, a=0.0nS, b=0.06nA, Vr = -55mV)
            return E1
        elseif SNN.AdExParameter(;El=-70mV, τm = 200ms, τe = 100ms, a=0.0nS, b=0.005nA, Vr = -55mV)
            return E2
        else  SNN.AdExParameter(;El=-70mV, τm = 5.0ms, τe = 100ms, a=-0.5nS, b=0.007nA, Vr = -46mV)
            return E3
        end
    end
end

#E.param

#E.param.R

#E.param.b

# Monitor neuron's voltage and adaptation variable
##
"Tonic firing: neuron fires at a const. rate in response to a steady input current. 
This causes the AdEx neuron to generate a regular series of action potentials that are separated by constant intervals."
# Simulate the neuron with a low amplitude step current
#param1=SNN.AdExParameter(El=-70mV, τm = 20ms, τe = 30ms, a=0.0nS, b=0.06nA, Vr = -55mV)
E1= initialize_pop()
E.I=[19]
SNN.sim!([E1], []; duration = 300ms)
p1 = plot(SNN.vecplot(E, :w),SNN.vecplot(E, :v), label = "Tonic firing")
##

"Adapting firing: neuron initially generates a high-frequency burst of action potentials in response to a current stimulus, but the firing rate then slows down and becomes more irregular over time. 
This is observed as a progressive decrease in the instantaneous firing rate in response to a const. current input."
# Simulate neuron with a high amplitude step current
#param2=SNN.AdExParameter2(El=-70mV, τm = 200ms, τe = 100ms, a=0.0nS, b=0.005nA, Vr = -55mV)
E2= initialize_pop()
E.I=[75]
SNN.sim!([E2], []; duration = 300ms)
p2 = plot(SNN.vecplot(E, :w),SNN.vecplot(E, :v), label = "Adapting firing")
##

"Bursting firing: neuron generates a series of high-frequency action potentials that are separated by 'resting' periods. 
This is observed as a series of spikes followed by a period of silence/low-level activity. 
Bursting firing can be classified into different subtypes based on the number and duration of the bursts."

# Simulate neuron with a 'burst' of high amplitude current
#param3=SNN.AdExParameter(El=-70mV, τm = 5.0ms, τe = 100ms, a=-0.5nS, b=0.007nA, Vr = -46mV)
E3 = initialize_pop()
E.I = [50, -50, -50, 50, -50, 50, -50, 50, -50, 50, -50]
SNN.sim!([E3], []; duration = 300ms)
p3=plot(SNN.vecplot(E, :w),SNN.vecplot(E, :v), label = "Bursting firing")
##
plot(p1,p2,p3, layout=(3,1), size=(900,1200))

xlabel!("Time")

ylabel!("Membrane Potential (mV)")
plot!(title="Multiple Firing Patterns\n in AdEx Neuron Model")





using Interact
using Plots

# Define the adaptation current parameter
AC = 10pA

# Create a function to plot the neuron response
function plot_neuron_response(AC)
  # Create an AdEx neuron with the specified adaptation current
  E = SNN.AdEx(; param = SNN.AdExParameter(; El = -70mV, τm = 20ms, Vt = -50mV, Vr = -70mV, AC = AC))

  # Set the input current
  E.I = 100pA

  # Initialize the membrane potential
  E.v .= -60mV

  # Monitor the membrane potential
  SNN.monitor(E, [:v])

  # Simulate the neuron for 1000ms
  SNN.sim!([E], []; duration = 1000ms)

  # Get the time vector and membrane potential vector
  t = SNN.time(E)
  V = E.records[:v]

  # Plot the neuron response
  plot(t, V, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)", title = "Neuron Response with AC = $AC", legend = :bottomright, label = "V")
end

# Create an interactive plot
@interact
function plot_interactive(AC = 10pA)
  plot_neuron_response(AC)
end
t = SNN.time(E)
plot_neuron_response(t, V)


using Revise
using DrWatson
@quickactivate "SNNExamples"
using Base
using SpikingNeuralNetworks
using Plots

default(
  size=(600, 400),
  tickfontsize=12,
  guidefontsize=14,
  margin=10Plots.mm,
  titlefontsize=18,
)

SNN.@load_units
const pA = ampere / 1e12

function initialize_pop(firing_pattern)
    param = param_configurations[firing_pattern]
    E = SNN.AdEx(;N=1, param=param)
    SNN.monitor(E, [:v, :fire, :w])
    return E
end


# Define neuron firing patterns and their parameter configurations
firing_patterns = ["tonic", "adapting", "bursting"]
param_configurations = Dict(
    "tonic" => dict(El=-70mV, τm=20ms, τe=30ms, a=0.0nS, b=60pA, Vr=-55mV, Vt=-50mV, ΔT=2mV),
    "adapting" => dict(El=-70mV, τm=20ms, τe=100ms, a=0.0nS, b=50nA, Vr=-55mV, Vt=-50mV, ΔT=2mV),
    "bursting" => dict(El=-70mV, τm=5.0ms, τe=100ms, a=-0.5nS, b=70nA, Vr=-46, Vt=-50mV, ΔT=2mV)
)

# Run simulations and generate plots for each firing pattern
for firing_pattern in firing_patterns
    # Initialize neuron population
    E = initialize_pop(firing_pattern)

    # Set input current
    input_current = eval(firing_pattern == "tonic" ? 40 : 65)
    E.I = [input_current]

    # Set initial membrane potential
    E.v .= -70mV

    # Simulate the network
    SNN.sim!([E], []; duration = 200ms)

    # Plot adaptation current and membrane potential
    plot(SNN.vecplot(E, :w), SNN.vecplot(E, :v), label=firing_pattern, legend=false, layout=(2,1))

    # Log AdEx neuron parameters
    println("AdEx Neuron Parameters for ", firing_pattern, ":")
    for x in fieldnames(SNN.AdExParameter)
        println(x, "  ", getfield(E.param, x))
    end
end

# Combine plots into a single layout
plot(p1, p2, p3, layout=(3, 1), legend=true, size=(900, 1600))
plot!(; legend=:topright)
xlabel!("Time (ms)")
ylabel!("Membrane Potential (mV)")
plot!(title="Multiple Firing Patterns in AdEx Neuron Model")

function initialize_pop(firing_pattern)
    param = param_configurations[firing_pattern]
    E = SNN.AdEx(;N=1, param=param)
    SNN.monitor(E, [:v, :fire, :w])
    return E
end


#

using Revise
using DrWatson
@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots
default(
    size=(900, 600),
    tickfontsize=10,
    guidefontsize=10,
    margin=8Plots.mm,
    titlefontsize=15,
    titlefontcolor=:darklavender
)

SNN.@load_units

τm = 30ms   #membrane capacitance
Vt = -50mV  #initial threshold 
Vr = -70mV #reference resting potential
El = Vr  #resting membrane potential

Is = 10:5:25
Is  

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = τm, Vt = Vt, Vr = Vr))

SNN.monitor(E, [:v,:fire])
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
#[plot!(x->Is[i],2000:6000) for i in eachindex(Is)]
plot!(x->0,6000:10000)
plot!(q)

plot(1:10:4000,my_Is[:,:]', label=string(Is))

f= (plot(plot(SNN.vecplot(E,:v),
q, layout=(2,1))))


using Revise
using DrWatson
@quickactivate "SNNExamples"

using SpikingNeuralNetworks
using Plots
default(
  size=(900, 600),
  tickfontsize=10,
  guidefontsize=10,
  margin=8Plots.mm,
  titlefontsize=15,
  titlefontcolor=:darklavender
)

SNN.@load_units

τm = 30ms  #membrane capacitance
Vt = -50mV  #initial threshold
Vr = -70mV #reference resting potential
El = Vr  #resting membrane potential

Is = 10:5:25
Is

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = τm, Vt = Vt, Vr = Vr))

SNN.monitor(E, [:v,:fire])
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

[plot!(x->Is[i],2000:6000, label=string(Is[i]) , legend = true, xlabel = "Time (ms)" ) for i in eachindex(Is)]

plot!(x->0,6000:10000)

f= plot!(
plot(SNN.vecplot(E,:v),
xlabel = "Time (ms)",
ylabel = "Membrane Potential (mV)",
title = "AdEx neuron in response to stochastic (noisy) input currents"),
q, layout=(2,1)
)
