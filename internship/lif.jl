#Simulation of the Leaky Integrate-and-Fire neural model
"[Integrate-And-Fire Neuron] https://neuronaldynamics.epfl.ch/online/Ch1.S3.html" 
"https://celltypes.brain-map.org/experiment/electrophysiology/474626527"

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

#Define neuron model parameters
τm = 30ms   #membrane capacitance
Vt = -50mV  #initial threshold 
Vr = -70mV #reference resting potential
El = Vr  #resting membrane potential

Is = 10:5:25 #current input range in pA
Is  

##
#Define AdEx neuron

E = SNN.AdEx(; N = length(Is), param = SNN.AdExParameter(; El = -70mV, τm = τm, Vt = Vt, Vr = Vr))
E.I = Is
E.v .= -70mV #initial membrane potential is set to resting potential 
SNN.monitor(E, [:v,:fire])
SNN.sim!([E], []; duration = 400ms) #simulate spiking behavior

SNN.vecplot(E, :v)
xlabel!("Time (ms)")
ylabel!("Membrane Potential (mV)")
plot!(title = "AdEx Neuron Model")


