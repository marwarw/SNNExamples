#Multiple firing patterns in the AdEx neuron model
"https://neuronaldynamics.epfl.ch/online/Ch6.S1.html"

using Revise
using DrWatson
@quickactivate "SNNExamples"

import SpikingNeuralNetworks
using Plots

default(
    size=(900, 600),
    tickfontsize=13,
    guidefontsize=18,
    margin=10Plots.mm,
    titlefontsize=18,
)

SNN.@load_units

const pA = ampere / 1e12

#AdEx neuron
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
    SNN.monitor(E, [:v, :fire, :w])
    return E
end


#Tonic firing:
"neuron fires at a const. rate in response to steady input current. 
=> AdEx neuron generates series of action potentials separated by constant intervals."

E1 = initialize_pop("tonic")
E1.I = [40]
SNN.sim!([E1], []; duration = 500ms)
p1 = plot(SNN.vecplot(E1, :w), SNN.vecplot(E1, :v), label="Tonic firing", legend=false)

#Adapting firing:
"neuron initially generates a high-frequency burst of action potentials 
in response to a current stimulus, then firing rate slows down and becomes more irregular over time. 
=> observed as a progressive decrease in the instantaneous firing rate in response to const. current input."

E2 = initialize_pop("adapting")
E2.I = [65]
SNN.sim!([E2], []; duration = 300ms)
p2 = plot(SNN.vecplot(E2, :w), SNN.vecplot(E2, :v), label="Adapting firing", legend=false)

#Bursting firing:
"neuron generates a series of high-frequency action potentials separated by resting periods. 
=> observed as a series of spikes followed by a period of silence/low-level activity. 
=> can be classified into different subtypes based on the number and duration of the bursts."

E3 = initialize_pop("bursting")
E3.I = [65,-65,65,-65,65,-65,65,-65]
SNN.sim!([E3], []; duration = 30ms)
p3 = plot(SNN.vecplot(E3, :w), SNN.vecplot(E3, :v), label="Bursting firing", legend=false)


plot(p1, p2, p3, layout=(3, 1), legend=true, size=(900, 1200))
plot!(; legend=:topright)
xlabel!("Time (ms)")
ylabel!("Membrane Potential (mV)")
plot!(title="Multiple Firing Patterns \n in AdEx Neuron Model")








