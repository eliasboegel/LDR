import Pkg
Pkg.add("DataFrames")
Pkg.add("Parameters")

include("../LDR.jl")
using .LDR

configurations::Vector{Params} = [
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=220e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=210e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=200e3, ablation_time=90)
    #
    # REFINEMENT LEVEL 1
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=300e3, ablation_time=87.5)
    #
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=77.5)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=87.5)
    Params(t0=7 * 24 * 3600, range=295e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=77.5)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=290e3, ablation_time=87.5)
    #
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=77.5)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=87.5)
    Params(t0=7 * 24 * 3600, range=285e3, ablation_time=90)
    #
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=77.5)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=280e3, ablation_time=87.5)
    #
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=77.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=80)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=82.5)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=85)
    Params(t0=7 * 24 * 3600, range=275e3, ablation_time=87.5)
    #
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=270e3, ablation_time=77.5)
    #
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=57.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=70)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=75)
    Params(t0=7 * 24 * 3600, range=265e3, ablation_time=77.5)
    #
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=55)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=57.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=72.5)
    Params(t0=7 * 24 * 3600, range=260e3, ablation_time=75)
    #
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=52.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=55)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=57.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=67.5)
    Params(t0=7 * 24 * 3600, range=255e3, ablation_time=70)
    #
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=47.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=52.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=55)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=57.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=65)
    Params(t0=7 * 24 * 3600, range=250e3, ablation_time=67.5)
    #
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=32.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=35)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=37.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=42.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=45)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=47.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=52.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=55)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=57.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=60)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=62.5)
    Params(t0=7 * 24 * 3600, range=245e3, ablation_time=65)
    #
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=32.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=35)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=37.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=42.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=45)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=47.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=52.5)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=55)
    Params(t0=7 * 24 * 3600, range=240e3, ablation_time=57.5)
    #
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=20)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=32.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=35)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=37.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=42.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=45)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=47.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=50)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=52.5)
    Params(t0=7 * 24 * 3600, range=235e3, ablation_time=55)
    #
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=22.5)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=25)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=32.5)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=35)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=37.5)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=42.5)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=45)
    Params(t0=7 * 24 * 3600, range=230e3, ablation_time=47.5)
    #
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=27.5)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=30)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=32.5)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=35)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=37.5)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=40)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=42.5)
    Params(t0=7 * 24 * 3600, range=225e3, ablation_time=45)
]

run_configurations(configurations, "runs.csv")