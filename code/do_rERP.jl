include("rERP.jl");

elec = [:Fp1, :Fp2, :F7, :F3, :Fz, :F4, :F8, :FC5, :FC1, :FC2, :FC6, :C3,
        :Cz,  :C4,  :CP5, :CP1, :CP2, :CP6, :P7,  :P3,  :Pz,  :P4,  :P8,  :O1,  :Oz,  :O2];

# Within-subjects regression
models = make_models([:Subject, :Timestamp], [:Item, :Condition], elec, [:Intercept, :Plaus, :Cloze_distractor]);
dt = process_data("../data/adbc23_erp.csv", false, models, invert_preds=[:Plaus]);

@time fit_models(dt, models, "rERPs_Plaus_Clozedist");

# Across-subjects regression, yielding a single-tvalue for each electrode and time-step
dts = dt;
dts.Subject = ones(nrow(dts));
@time fit_models(dts, models, "rERPs_Plaus_Clozedist_across");

# Predict EEG by ReadingTimes (from separate Exp)
models = make_models([:Subject, :Timestamp], [:Item, :Condition], elec, [:Intercept, :ReadingTime]);
@time dtrt = process_data("../data/adbc23_erp.csv", false, models);

@time fit_models(dtrt, models, "rERPS_ReadingTime");