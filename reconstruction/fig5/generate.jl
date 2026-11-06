include(joinpath(@__DIR__, "..", "study_io.jl"))

paper_settings = (
    mu=3.0, nmax=150, T=5.0, dt=0.01, tau=0.05, eta=0.9, b=0.0,
    lambda=1e-6, solver_tol=5e-13, seed=20261016, n_runs=30, R_fixed=10_000,
    R_values=sort(unique(vcat(10_000, round.(Int, 10.0 .^ range(2, 5; length=9))))),
    eta_values=[0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.0],
    mean_values=[1.0, 2.0, 3.0, 5.0, 10.0, 20.0], tail_tolerance=1e-10,
    families=("coherent", "thermal"),
)

# 1. Main paper study: kernels and observations, with no optimization yet.
# R subsets and eta thinning retain precisely the original coupling and seeds.
function paper_observations(s)
    check_detector_settings(s)
    (; mu, nmax, T, dt, tau, eta, b, seed, n_runs, R_fixed,
       R_values, eta_values, mean_values, tail_tolerance, families) = s
    @assert issorted(R_values) && allunique(R_values) && all(R_values .> 0)
    @assert R_fixed in R_values && all(0 .< eta_values .<= 1)
    @assert all(mean_values .> 0) && all(f in ("coherent", "thermal") for f in families)
    edges = collect(range(0.0, T; length=round(Int, T / dt) + 1))
    B = expected_background(edges, b, tau)
    all_means = sort(unique(vcat(mu, mean_values)))
    cutoffs = Dict(occupation => max(nmax,
        Int(quantile(Poisson(occupation), 1-tail_tolerance)),
        Int(quantile(Geometric(1/(1+occupation)), 1-tail_tolerance)))
        for occupation in all_means)
    println("Building paper kernels through n = $(maximum(values(cutoffs)))...")
    _, K_bank = intensity_kernels(maximum(values(cutoffs)), edges, tau)
    cases = Dict{Tuple{String,Float64},Any}()
    for occupation in all_means, family in families
        cutoff = cutoffs[occupation]
        law = family == "coherent" ? Poisson(occupation) : Geometric(1/(1+occupation))
        p_true = pdf.(law, 0:cutoff)
        p_true ./= sum(p_true)
        truth = number_moments(p_true)
        true_variance = family == "coherent" ? occupation : occupation*(1+occupation)
        @assert ccdf(law, cutoff) <= 1.01tail_tolerance
        @assert isapprox(truth.mean, occupation; rtol=1e-6)
        @assert isapprox(truth.variance, true_variance; rtol=1e-6)
        @assert isapprox(truth.g2, family == "coherent" ? 1.0 : 2.0; atol=1e-6)
        cases[(family, occupation)] = (; n=collect(0:cutoff), p_true, truth)
    end

    n_bins = length(edges)-1
    R_index = findfirst(==(R_fixed), R_values)
    R_data = Dict(f => zeros(n_bins, length(R_values), n_runs) for f in families)
    eta_data = Dict(f => zeros(n_bins, length(eta_values), n_runs) for f in families)
    N_data = Dict(f => zeros(n_bins, length(mean_values), n_runs) for f in families)
    jobs = [(family_index, family, run) for (family_index, family) in enumerate(families)
            for run in 1:n_runs]
    Threads.@threads for job in eachindex(jobs)
        family_index, family, run = jobs[job]
        rng = MersenneTwister(seed + 100_000run + 1_000family_index)
        pool = photon_record(rng, cases[(family, mu)].p_true, maximum(R_values),
            1.0, 0.0, T; keep_repetitions=true)
        marks = rand(rng, length(pool.events))
        background = photon_record(rng, [1.0], maximum(R_values),
            1.0, b, T; keep_repetitions=true)
        total = zeros(n_bins)
        previous_end = previous_background_end = 0
        for (r_index, repetitions) in enumerate(R_values)
            last_end = pool.preparation_ends[repetitions]
            indices = previous_end+1:last_end
            selected = view(pool.events, indices)[view(marks, indices) .< eta]
            total .+= detector_bins_sorted(selected, edges, tau)
            background_end = background.preparation_ends[repetitions]
            total .+= detector_bins_sorted(view(background.events,
                previous_background_end+1:background_end), edges, tau)
            R_data[family][:, r_index, run] = total / repetitions
            previous_end = last_end
            previous_background_end = background_end
        end
        fixed_end = pool.preparation_ends[R_fixed]
        fixed_signal = view(pool.events, 1:fixed_end)
        fixed_marks = view(marks, 1:fixed_end)
        background_y = detector_bins_sorted(view(background.events,
            1:background.preparation_ends[R_fixed]), edges, tau) / R_fixed
        for (i, efficiency) in enumerate(eta_values)
            eta_data[family][:, i, run] = efficiency == eta ? R_data[family][:, R_index, run] :
                detector_bins_sorted(fixed_signal[fixed_marks .< efficiency], edges, tau) / R_fixed + background_y
        end
        for (i, occupation) in enumerate(mean_values)
            if occupation == mu
                N_data[family][:, i, run] = R_data[family][:, R_index, run]
            else
                rng_N = MersenneTwister(seed + 100_000run + 1_000family_index + i)
                events = photon_record(rng_N, cases[(family, occupation)].p_true, R_fixed, eta, b, T)
                N_data[family][:, i, run] = detector_bins_sorted(events, edges, tau) / R_fixed
            end
        end
        println("Generated paper data: $family, run $run/$n_runs")
    end
    return (; settings=s, edges, B, K_bank, cases, R_data, eta_data, N_data)
end

# 2. Fit cached paper observations. No photons or kernels are regenerated here.
function fit_paper(observations; lambda=observations.settings.lambda,
                   solver_tol=observations.settings.solver_tol)
    s = merge(observations.settings, (; lambda, solver_tol))
    (; mu, eta, n_runs, R_fixed, R_values, eta_values, mean_values, families) = s
    (; B, K_bank, cases, R_data, eta_data, N_data) = observations
    R_index = findfirst(==(R_fixed), R_values)
    R_results = Dict(f => Matrix{Any}(undef, length(R_values), n_runs) for f in families)
    eta_results = Dict(f => Matrix{Any}(undef, length(eta_values), n_runs) for f in families)
    N_results = Dict(f => Matrix{Any}(undef, length(mean_values), n_runs) for f in families)
    jobs = [(f, run) for f in families for run in 1:n_runs]
    Threads.@threads for job in eachindex(jobs)
        family, run = jobs[job]
        K = view(K_bank, :, 1:length(cases[(family, mu)].p_true))
        for i in eachindex(R_values)
            fit = reconstruct(K, R_data[family][:, i, run], B, eta; lambda, solver_tol)
            R_results[family][i, run] = (; fit, moments=number_moments(fit.p))
        end
        for (i, efficiency) in enumerate(eta_values)
            if efficiency == eta
                eta_results[family][i, run] = R_results[family][R_index, run]
            else
                fit = reconstruct(K, eta_data[family][:, i, run], B, efficiency; lambda, solver_tol)
                eta_results[family][i, run] = (; fit, moments=number_moments(fit.p))
            end
        end
        for (i, occupation) in enumerate(mean_values)
            if occupation == mu
                N_results[family][i, run] = R_results[family][R_index, run]
            else
                K_N = view(K_bank, :, 1:length(cases[(family, occupation)].p_true))
                fit = reconstruct(K_N, N_data[family][:, i, run], B, eta; lambda, solver_tol)
                N_results[family][i, run] = (; fit, moments=number_moments(fit.p))
            end
        end
        println("Fitted paper data: $family, run $run/$n_runs")
    end
    clean_keys = [(f, occupation, efficiency) for f in families
        for occupation in sort(unique(vcat(mu, mean_values)))
        for efficiency in (occupation == mu ? sort(unique(vcat(eta, eta_values))) : [eta])]
    clean_fits = Vector{Any}(undef, length(clean_keys))
    Threads.@threads for i in eachindex(clean_keys)
        family, occupation, efficiency = clean_keys[i]
        case = cases[(family, occupation)]
        K = view(K_bank, :, 1:length(case.p_true))
        fit = reconstruct(K, efficiency*K*case.p_true+B, B, efficiency; lambda, solver_tol)
        clean_fits[i] = (; fit, moments=number_moments(fit.p))
    end
    clean_results = Dict(zip(clean_keys, clean_fits))
    summary_results = Dict{Tuple{Symbol,String,Symbol},Any}()
    for (scan, collection) in ((:R, R_results), (:eta, eta_results), (:N, N_results)),
        family in families, quantity in (:mean, :variance, :g2)
        fits = collection[family]
        samples = [getproperty(fits[i, run].moments, quantity)
                   for i in axes(fits, 1), run in axes(fits, 2)]
        @assert all(isfinite, samples)
        center = vec(mean(samples; dims=2))
        spread = vec(std(samples; dims=2, corrected=true))
        truth = [getproperty(cases[(family, scan == :N ? mean_values[i] : mu)].truth,
                    quantity) for i in axes(samples, 1)]
        rmse = sqrt.(vec(mean((samples .- truth).^2; dims=2)))
        summary_results[(scan, family, quantity)] = (;
            samples, center, spread, lower=center-spread, upper=center+spread, truth, rmse)
    end
    distribution_results = Dict{String,Any}()
    for family in families
        p = hcat([R_results[family][R_index, run].fit.p for run in 1:n_runs]...)
        average, spread = vec(mean(p; dims=2)), vec(std(p; dims=2, corrected=true))
        @assert isapprox(sum(average), 1; atol=1e-8)
        distribution_results[family] = (; average, spread, lower=average-spread, upper=average+spread)
    end
    return (; settings=s, R_results, eta_results, N_results, clean_results,
        summary_results, distribution_results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia --project=. fig5/generate.jl [output.jld2]")
    output = isempty(ARGS) ? joinpath(@__DIR__, "output", "fig5_results.jld2") : abspath(ARGS[1])
    generate_study("paper", paper_settings, paper_observations, fit_paper, output)
end
