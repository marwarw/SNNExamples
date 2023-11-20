#AEIF: Simulate adex model with spike after depolarization current and adaptive threshold
#R. Brette and W. Gerstner. Adaptive exponential integrate-and-fire model as an effective description of neuronal activity. J. Neurophysiol., 94:3637 – 3642, 2005.

function aEIF(u, w, w_tail, V_T, I, counter)
    # Model parameters
    th = 20            # [mV] spike threshold
    C = 281            # [pF] membrane capacitance
    g_L = 30           # [nS] membrane conductance
    E_L = -70.6        # [mV] resting voltage
    VT_rest = -50.4    # [mV] resetting voltage
    Delta_T = 2        # [mV] exponential parameters
    tau_w = 144        # [ms] time constant for adaptation variable w
    a = 4              # [nS] adaptation coupling constant
    b = 0.0805         # [nA] spike-triggered adaptation
    w_jump = 400       # spike after depolarization
    tau_wtail = 40     # [ms] time constant for spike after depolarization
    tau_VT = 50        # [ms] time constant for VT
    VT_jump = 20       # adaptive threshold

    if counter == 2
        u = E_L + 15 + 6.0984
        w += b
        w_tail = w_jump
        counter = 0
        V_T = VT_jump + VT_rest
    end

    udot = 1 / C * (-g_L * (u - E_L) + g_L * Delta_T * exp((u - V_T) / Delta_T) - w + w_tail + I)
    wdot = 1 / tau_w * (a * (u - E_L) - w)
    u += udot
    w += wdot
    w_tail -= w_tail / tau_wtail
    V_T = VT_rest / tau_VT + (1 - 1 / tau_VT) * V_T

    if counter == 1
        counter = 2
        u = 29.4 + 3.462
        w -= wdot
    end

    if u > th && counter == 0
        u = 29.4
        counter = 1
    end

    return u, w, w_tail, V_T, counter
end



using Plots

function simulate_plasticity(par, rho, DeltaTau)
    # Simulation of the model
    nFreq = length(rho)
    w = zeros(Float64, nFreq)
    for i in 1:nFreq
        w[i] = (VoTri(rho[i], DeltaTau, par) - 0.5) / 0.5
    end
    return w
end

# Parameters
par = [0.00014, 0.00008, 7, 15, 10]  # parameters of the plasticity model
rho = [0.1, 10, 20, 40, 50]          # frequency of the pairing for the simulation
DeltaTau = 10                        # time lag

# Pre-post pairing
w_prepost = simulate_plasticity(par, rho, DeltaTau)

# Plot pre-post pairing
plot(rho, w_prepost * 100 + 100, linewidth = 3, label = "Pre-Post Pairing")
xlabel!("ρ [Hz]")
ylabel!("Normalized Weight [%]")
title!("Voltage-Triplet Plasticity (Pre-Post Pairing)")


# Post-pre pairing
w_postpre = simulate_plasticity(par, rho, -DeltaTau)

# Plot post-pre pairing
plot!(rho, w_postpre * 100 + 100, linewidth = 2, color = :red, label = "Post-Pre Pairing")
