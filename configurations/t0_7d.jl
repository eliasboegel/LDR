include("../LDR.jl")
using .LDR

configurations::Vector{Params} = [
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=25)
    #
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=25)
    #
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=25)
    #
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=25)
    #
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=25)
    #
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=90)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=25)
]

run_configurations(configurations, "runs.csv")