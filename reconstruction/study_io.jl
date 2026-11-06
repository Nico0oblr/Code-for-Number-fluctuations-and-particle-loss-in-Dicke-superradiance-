include(joinpath(@__DIR__, "reconstruction.jl"))
using Statistics

# Same exponential detector response as detector_bins, but each sorted click is
# visited once. With tau=0 this reduces to ordinary counts in [left, right).
function detector_bins_sorted(events, edges, tau)
    clicks = sort(events)
    y = zeros(length(edges) - 1)
    tail = 0.0
    event_index = 1
    for j in eachindex(y)
        left, right = edges[j], edges[j+1]
        decay = tau == 0 ? 0.0 : exp(-(right-left) / tau)
        y[j] = tail * (1-decay)
        tail *= decay
        while event_index <= length(clicks) && clicks[event_index] < right
            s = clicks[event_index]
            if tau == 0
                y[j] += 1
            else
                y[j] += -expm1(-(right-s) / tau)
                tail += exp(-(right-s) / tau)
            end
            event_index += 1
        end
    end
    return y
end

function check_detector_settings(s)
    @assert s.mu > 0 && s.nmax >= 1 && s.n_runs >= 2
    @assert s.T > 0 && s.dt > 0 && s.tau >= 0 && s.b >= 0
    @assert 0 < s.eta <= 1 && isapprox(round(s.T / s.dt) * s.dt, s.T)
end

# Save binned observations first so a completed simulation can be refit.
function generate_study(study, settings, observe, fit, output)
    ispath(output) && error("$output already exists; choose another output filename")
    mkpath(dirname(output))
    BLAS.set_num_threads(1)
    println("Study: $study; Julia threads: $(Threads.nthreads())")
    provenance = (; study, julia_version=string(VERSION), threads=Threads.nthreads(), schema_version=1)
    observations = observe(settings)
    jldsave(output; observations, provenance)
    println("Saved observations to $output; starting fits...")
    fits = fit(observations)
    jldopen(output, "a") do file
        file["fits"] = fits
    end
    println("Saved $output")
end
