# A small, independent reference calculation. Read this file from top to bottom.
# Load the function definitions from your notebook with include("reconstruction.jl").
#
# Units: Gamma = 1. p[n+1] means P(N=n). K and y contain BIN INTEGRALS,
# not rates: y[j]/dt would be a rate for a uniform bin of width dt.
# This script imports none of the old reconstruction code.

using OrdinaryDiffEq, Distributions, Random, LinearAlgebra, JLD2
using JuMP, Clarabel

# 1. Forward model: solve the population ODE once for each fixed particle number.
# P_m' = lambda_(m+1)*P_(m+1) - lambda_m*P_m,
# lambda_m = m*(n-m+1), P_n(0)=1, I_n = sum_m lambda_m*P_m.
# The detector satisfies tau*z' = I_n-z and C'=z, initially z=C=0.
# Thus K[j,n+1] = C(t[j+1])-C(t[j]) includes BOTH tau and dt.
function intensity_kernels(nmax, edges, tau)
    I = zeros(length(edges), nmax + 1)
    K = zeros(length(edges) - 1, nmax + 1)
    for n in 1:nmax
        z_index, C_index = n + 2, n + 3
        u0 = zeros(n + 3)
        u0[n + 1] = 1.0
        function rhs!(du, u, parameters, t)
            fill!(du, 0.0)
            intensity = 0.0
            for m in 1:n
                flux = m * (n - m + 1) * u[m + 1]
                du[m + 1] -= flux
                du[m] += flux
                intensity += flux
            end
            if tau > 0
                du[z_index] = (intensity - u[z_index]) / tau
                du[C_index] = u[z_index]
            else
                du[C_index] = intensity
            end
        end
        solution = solve(ODEProblem(rhs!, u0, (0.0, edges[end])), Tsit5();
                         saveat=edges, save_everystep=false,
                         reltol=1e-10, abstol=1e-12)
        string(solution.retcode) == "Success" || error("Population ODE failed")
        for j in eachindex(edges)
            I[j, n + 1] = sum(m * (n-m+1) * solution.u[j][m+1] for m in 1:n)
        end
        K[:, n + 1] = diff([u[C_index] for u in solution.u])
    end
    return I, K
end

# 2. Noisy photon records: draw N, then classical exponential waiting times.
# Keep each emitted photon with probability eta. Add Poisson(b*T) background
# clicks, uniformly distributed in time. b is the REGISTERED background rate.
# These random events already produce counting noise; add no further noise.
function photon_record(rng, p, R, eta, b, T; keep_repetitions=false)
    events = Float64[]
    preparation_ends = keep_repetitions ? zeros(Int, R) : Int[]
    for repetition in 1:R
        n = rand(rng, Categorical(p)) - 1
        t = 0.0
        for m in n:-1:1
            t += randexp(rng) / (m * (n - m + 1))
            t >= T && break
            rand(rng) < eta && push!(events, t)
        end
        for background_event in 1:rand(rng, Poisson(b * T))
            push!(events, rand(rng) * T)
        end
        keep_repetitions && (preparation_ends[repetition] = length(events))
    end
    # The first r preparations occupy events[1:preparation_ends[r]].
    return keep_repetitions ? (; events, preparation_ends) : events
end

# 3. Detector and bins: each click at s contributes the exact area of
# h(t-s)=exp(-(t-s)/tau)/tau for t>=s in [left,right), right-left=dt.
# This explicit double loop prioritizes readability over large-scan speed.
function detector_bins(events, edges, tau)
    y = zeros(length(edges) - 1)
    for s in events, j in eachindex(y)
        left, right = edges[j], edges[j + 1]
        if tau == 0
            y[j] += left <= s < right
        elseif s < left
            y[j] += exp(-(left-s)/tau) * (-expm1(-(right-left)/tau))
        elseif s < right
            y[j] += -expm1(-(right-s)/tau)
        end
    end
    return y
end

# Expected background BIN AREA per preparation, with a zero initial filter:
# B_j = b*[dt - tau*(exp(-left/tau)-exp(-right/tau))].
function expected_background(edges, b, tau)
    return [tau == 0 ? b * (right-left) :
            b * ((right-left) - tau * (exp(-left/tau)-exp(-right/tau)))
            for (left, right) in zip(edges[1:end-1], edges[2:end])]
end

# 4. Reconstruction: library optimization, with the objective and constraints
# written explicitly. Residual variables avoid explicitly forming K'*K.
# A single common scale improves numerical units without reweighting bins.
# Minimize ||eta*K*p+B-y||^2 + lambda*sum((p[n+1]-2p[n]+p[n-1])^2).
# lambda=0 reproduces the unregularized fit. Positive lambda prefers smooth p.
# solver_tol controls objective-gap accuracy; feasibility tolerance is 10*solver_tol.
# Tighter requests also refine the internal linear solves, without changing the objective.
function reconstruct(K, y, B, eta; lambda::Real=0.0, solver_tol::Real=1e-12)
    isfinite(lambda) && lambda >= 0 || error("lambda must be finite and nonnegative")
    isfinite(solver_tol) && solver_tol > 0 || error("solver_tol must be finite and positive")
    J, L = size(K)
    scale = max(norm(y - B), eps(Float64))
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol_feas", 10solver_tol)
    set_optimizer_attribute(model, "tol_gap_abs", solver_tol)
    set_optimizer_attribute(model, "tol_gap_rel", solver_tol)
    numerical_regularization = solver_tol < 1e-12 ? 100solver_tol : 1e-8
    if solver_tol < 1e-12
        set_optimizer_attribute(model, "max_iter", 400)
        set_optimizer_attribute(model, "static_regularization_constant", numerical_regularization)
        set_optimizer_attribute(model, "iterative_refinement_abstol", solver_tol / 10)
        set_optimizer_attribute(model, "iterative_refinement_reltol", solver_tol / 10)
        set_optimizer_attribute(model, "iterative_refinement_max_iter", 30)
        set_optimizer_attribute(model, "dynamic_regularization_eps", solver_tol / 100)
        set_optimizer_attribute(model, "dynamic_regularization_delta", sqrt(solver_tol))
    end
    @variable(model, p[1:L] >= 0)
    @variable(model, residual[1:J])
    @constraint(model, sum(p) == 1)
    @constraint(model, [j in 1:J],
        scale * residual[j] == eta * sum(K[j,n] * p[n] for n in 1:L) + B[j] - y[j])
    # Divide the ENTIRE objective by scale^2, so lambda has the meaning above.
    @objective(model, Min,
        sum(residual[j]^2 for j in 1:J) +
        (lambda / scale^2) * sum((p[n+1] - 2p[n] + p[n-1])^2 for n in 2:(L-1)))
    optimize!(model)
    # A nearly singular KKT system can stall at the requested accuracy.
    # Retry with smaller numerical regularization, keeping the objective and tolerances.
    if termination_status(model) == JuMP.MOI.ALMOST_OPTIMAL
        set_optimizer_attribute(model, "static_regularization_constant", numerical_regularization / 100)
        optimize!(model)
    end
    termination_status(model) == JuMP.MOI.OPTIMAL ||
        error("Reconstruction failed at solver_tol=$solver_tol: $(termination_status(model))")
    p_fit = value.(p)
    return (p=p_fit, lambda=Float64(lambda), solver_tol=Float64(solver_tol),
            roughness=sum(abs2, diff(diff(p_fit))),
            status=string(termination_status(model)),
            residual_norm=norm(eta * K * p_fit + B - y),
            normalization_error=abs(sum(p_fit)-1), minimum_probability=minimum(p_fit))
end

# 5. Statistics of the INITIAL particle number, not optical coincidences.
function number_moments(p)
    n = collect(0:length(p)-1)
    mean_N = dot(n, p)
    variance_N = dot(n.^2, p) - mean_N^2
    g2_N = mean_N > 0 ? dot(n .* (n .- 1), p) / mean_N^2 : NaN
    return (mean=mean_N, variance=variance_N, g2=g2_N)
end
