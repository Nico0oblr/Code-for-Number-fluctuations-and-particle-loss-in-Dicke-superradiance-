include(joinpath(@__DIR__, "..", "study_io.jl"))

lambda_settings = (
    mu=5.0, nmax=150, T=20.0, dt=0.01, tau=0.05, eta=0.9, b=0.0,
    R=100_000, n_runs=30, seed=42, solver_tol=5e-13,
    lambda_values=vcat(0.0, 10.0 .^ range(-7, -1; length=25)),
    families=("coherent", "thermal", "fock", "even-cat"),
    moment_lambda=0.0, moment_orders=collect(2:4),
    distribution_lambdas=[0.0, 1e-5], distribution_n_display=40,
)

# 3. Smoothing appendix: four states with the same mean and cached noisy curves.
function lambda_observations(s)
    check_detector_settings(s)
    (; mu, nmax, T, dt, tau, eta, b, R, n_runs, seed, families) = s
    @assert R > 0 && isinteger(mu) && mu <= nmax
    @assert all(f in ("coherent", "thermal", "fock", "even-cat") for f in families)
    edges = collect(range(0.0, T; length=round(Int, T/dt)+1))
    n = collect(0:nmax)
    println("Building smoothing-study kernels through n = $nmax...")
    _, K = intensity_kernels(nmax, edges, tau)
    B = expected_background(edges, b, tau)
    lower, upper = mu, mu+1.0
    for _ in 1:80
        midpoint = (lower+upper)/2
        if midpoint*tanh(midpoint) < mu
            lower = midpoint
        else
            upper = midpoint
        end
    end
    cat_x = (lower+upper)/2
    records = Dict{String,Any}()
    for family in families
        if family == "fock"
            p_true = zeros(nmax+1)
            p_true[Int(mu)+1] = 1.0
        elseif family == "even-cat"
            law = Poisson(cat_x)
            even_probability = (1+exp(-2cat_x))/2
            @assert ccdf(law, nmax)/even_probability <= 1e-10
            p_true = [iseven(number) ? pdf(law, number)/even_probability : 0.0 for number in n]
        else
            law = family == "coherent" ? Poisson(mu) : Geometric(1/(1+mu))
            @assert ccdf(law, nmax) <= 1e-10 "Increase nmax for $family"
            p_true = pdf.(law, n)
        end
        p_true ./= sum(p_true)
        @assert isapprox(number_moments(p_true).mean, mu; rtol=1e-7)
        records[family] = (; p_true, truth=number_moments(p_true).g2,
            expected_y=eta*K*p_true+B, y=zeros(length(edges)-1, n_runs))
    end
    jobs = [(i, f, run) for (i, f) in enumerate(families) for run in 1:n_runs]
    Threads.@threads for job in eachindex(jobs)
        family_index, family, run = jobs[job]
        rng = MersenneTwister(seed+10_000family_index+run)
        events = photon_record(rng, records[family].p_true, R, eta, b, T)
        records[family].y[:, run] = detector_bins_sorted(events, edges, tau)/R
        println("Generated smoothing data: $family, run $run/$n_runs")
    end
    return (; settings=s, edges, n, K, B, records)
end

function number_correlations(p, orders)
    numbers = Float64.(0:length(p)-1)
    mean_N = dot(numbers, p)
    mean_N > 0 || error("Normalized correlations require a positive mean")
    return [sum(p .* [prod(number-j for j in 0:k-1) for number in numbers])/mean_N^k
            for k in orders]
end

# All p samples and g2/g3/g4 summaries are saved for EVERY lambda. Choosing a
# different displayed lambda in the plot notebook never calls the optimizer.
function fit_lambda(observations; lambda_values=observations.settings.lambda_values,
                    solver_tol=observations.settings.solver_tol)
    s = merge(observations.settings, (; lambda_values=copy(lambda_values), solver_tol))
    (; eta, n_runs, families, moment_orders) = s
    (; K, B, records) = observations
    @assert all(isfinite, lambda_values) && all(lambda_values .>= 0)
    results = Dict(f => Vector{Any}(undef, length(lambda_values)) for f in families)
    jobs = [(f, i) for f in families for i in eachindex(lambda_values)]
    Threads.@threads for job in eachindex(jobs)
        family, i = jobs[job]
        smoothing = lambda_values[i]
        record = records[family]
        probabilities = zeros(size(K, 2), n_runs)
        statuses = Vector{String}(undef, n_runs)
        residuals = zeros(n_runs)
        for run in 1:n_runs
            fit = reconstruct(K, record.y[:, run], B, eta; lambda=smoothing, solver_tol)
            probabilities[:, run] = fit.p
            statuses[run], residuals[run] = fit.status, fit.residual_norm
        end
        clean_fit = reconstruct(K, record.expected_y, B, eta; lambda=smoothing, solver_tol)
        samples = hcat([number_correlations(p, moment_orders) for p in eachcol(probabilities)]...)
        g2 = [number_moments(p).g2 for p in eachcol(probabilities)]
        @assert all(isfinite, samples)
        moments = (; samples, truth=number_correlations(record.p_true, moment_orders),
            average=vec(mean(samples; dims=2)), spread=vec(std(samples; dims=2, corrected=true)))
        distribution = (; lambda=smoothing, average=vec(mean(probabilities; dims=2)),
            spread=vec(std(probabilities; dims=2, corrected=true)))
        results[family][i] = (; lambda=smoothing, solver_tol, g2, p=probabilities,
            average=mean(g2), spread=std(g2; corrected=true),
            noiseless=number_moments(clean_fit.p).g2, clean_fit, moments, distribution,
            statuses, residuals)
        println("Fitted smoothing data: $family, lambda=$smoothing")
    end
    return (; settings=s, results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia --project=. figS4/generate.jl [output.jld2]")
    output = isempty(ARGS) ? joinpath(@__DIR__, "output", "figS4_results.jld2") : abspath(ARGS[1])
    generate_study("lambda", lambda_settings, lambda_observations, fit_lambda, output)
end
