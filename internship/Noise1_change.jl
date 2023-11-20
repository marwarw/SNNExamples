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

τm = 30ms
Vt = -50mV
Vr = -70mV
El = Vr
Is = 10:5:25 # in pA 
N_neurons = 5
E = SNN.AdEx(; N = N_neurons, param = SNN.AdExParameter(; El = El, τm = τm, Vt = Vt, Vr = Vr))
E.v .=El
in_rate = 200Hz 

# Combination 1: (σ1; μ1:μ2)
σ1 = 5
μ1 = 10
μ2 = 100
rfn1 = (n) -> rand(Normal(μ1, σ1), n)

in_rates1 = range(μ1, μ2, N_neurons) * Hz
P1, C1 = [], []
push!(P1, E)
for i in 1:N_neurons
    in_rate1 = in_rates1[i]
    D1 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate1))
    EE1 = SNN.SpikingSynapse(D1, E, :ge; rfn = rfn1, p = 1.0)
    EE1.W .= 0
    EE1.W[i] = σ1
    push!(C1, EE1)
    push!(P1, D1)
end

SNN.monitor(E, [:v, :fire])
SNN.sim!(P1, C1; duration = 1000ms)
v1 = SNN.getrecord(E, :v)
y1 = hcat(v1...)'

# Combination 2: (σ2; μ1:μ2)
σ2 = 50
rfn2 = (n) -> rand(Normal(μ1, σ2), n)
P2, C2 = [], []
push!(P2, E)
in_rates2 = range(μ1, μ2, N_neurons) * Hz

for i in 1:N_neurons
    in_rate2 = in_rates2[i]
    D2 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate2))
    EE2 = SNN.SpikingSynapse(D2, E, :ge; rfn = rfn2, p = 1.0)
    EE2.W .= 0
    EE2.W[i] = 10
    push!(C2, EE2)
    push!(P2, D2)
end

SNN.monitor(E, [:v, :fire])
SNN.sim!(P2, C2; duration = 1000ms)
v2 = SNN.getrecord(E, :v)
y2 = hcat(v2...)'

# Combination 3: (σ1:σ2; μ1)
rfn3 = (n,ξ) -> rand(ξ % 2 == 0 ? Normal(μ1, σ1) : Normal(μ1, σ2), n)

σ2
plot(histogram(rand(Normal(10,σ1),1000)),
    histogram(rand(Normal(10,σ2),1000)), layout=(2,1), link=:x)

n=1000
plot(histogram(rfn3(n,1)),
    histogram(rfn3(n,2)), layout=(2,1), link=:x)


rfn3(10)
P3, C3 = [], []
push!(P3, E)
in_rates3 = repeat([in_rate], N_neurons)

for i in 1:N_neurons
    in_rate3 = in_rates3[i]
    _rfn3(n) = rfn3(n, i)
    D3 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate3))
    EE3 = SNN.SpikingSynapse(D3, E, :ge; rfn = _rfn3, p = 1.0)
    EE3.W .= 0
    EE3.W[i] = 10
    push!(C3, EE3)
    push!(P3, D3)
end

SNN.monitor(E, [:v, :fire])
SNN.sim!(P3, C3; duration = 1000ms)
v3 = SNN.getrecord(E, :v)
y3 = hcat(v3...)'

# Combination 4: (σ1:σ2; μ2)
rfn4 = (n) -> rand(n % 2 == 0 ? Normal(μ2, σ1) : Normal(μ2, σ2), n)
P4, C4 = [], []
push!(P4, E)
in_rates4 = repeat([in_rate], N_neurons)

for i in 1:N_neurons
    in_rate4 = in_rates4[i]
    D4 = SNN.Poisson(; N = 1, param = SNN.PoissonParameter(; rate = in_rate4))
    EE4 = SNN.SpikingSynapse(D4, E, :ge; rfn = rfn4, p = 1.0)
    EE4.W .= 0
    EE4.W[i] = 10
    push!(C4, EE4)
    push!(P4, D4)
end

SNN.monitor(E, [:v, :fire])
SNN.sim!(P4, C4; duration = 1000ms)
v4 = SNN.getrecord(E, :v)
y4 = hcat(v4...)'

# Plotting
x = 1:length(v1)
p1 = plot(x, y1, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)",
    title = "Combination 1: (σ1; μ1:μ2)", leg = :none,
    xaxis = ("Time (ms)", extrema(x)),
    yaxis = ("Membrane Potential (mV)", extrema(y1)), palette = (:vanimo10), lw = 4)

p2 = plot(x, y2, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)",
    title = "Combination 2: (σ2; μ1:μ2)", leg = :none,
    xaxis = ("Time (ms)", extrema(x)),
    yaxis = ("Membrane Potential (mV)", extrema(y2)), palette = (:vanimo10), lw = 4)

p3 = plot(x, y3, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)",
    title = "Combination 3: (σ1:σ2; μ1)", leg = :none,
    xaxis = ("Time (ms)", extrema(x)),
    yaxis = ("Membrane Potential (mV)", extrema(y3)), palette = (:vanimo10), lw = 4)

p4 = plot(x, y4, xlabel = "Time (ms)", ylabel = "Membrane Potential (mV)",
    title = "Combination 4: (σ1:σ2; μ2)", leg = :none,
    xaxis = ("Time (ms)", extrema(x)),
    yaxis = ("Membrane Potential (mV)", extrema(y4)), palette = (:vanimo10), lw = 4)

plot(p1, p2, p3, p4, layout = (2, 2), size = (900, 600))

p1

#heatmap
σ_values = [σ1, σ2]
μ_values = [μ1, μ2]
m_potentials = [y1, y2, y3, y4]
heatmap(μ_values, σ_values, m_potentials,
    xlabel = "Standard Deviation (σ)",
    ylabel = "Mean (μ)",
    title = "AdEx Neuron Model: σ vs. μ",
    color = :vanimo10,
    c = :blues,
    aspect_ratio = :equal,
    colorbar_title = "Membrane Potential (mV)"
)





