#Noise 1 (Synaptic noise) in tonic firing, with high σ, high rate
using Revise
using DrWatson
using Logging

@quickactivate "SNNExamples"

using Plots
using Interact

plotlyjs()
@manipulate for var1=1:20
    plot(x->var1*x)
end