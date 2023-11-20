using Revise
using DrWatson
@quickactivate "SNNExamples"

using Plots
using SpikingNeuralNetworks
SNN.@load_units

# Define the parameters
Is = [10, 30, 50, 70, 90]  # Increasing current step values

firing_rate = []  # Firing rate array

# Define the Rate neuron
E = SNN.Rate(; N = length(Is))

# Simulate and collect firing rates

EE = SNN.RateSynapse(E, E, σ = 5, p = 1.0)
    
SNN.monitor(E, [(:r)])

SNN.sim!([E], [EE]; duration = 1000ms)


# Plot the firing rates
#plot(Is, firing_rate, marker = :circle, xlabel = "Current Step (Is)", ylabel = "Firing Rate (spikes/s)", legend = false)

SNN.vecplot(E, :r)