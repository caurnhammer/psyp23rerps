include("lmerRT.jl");

## Data Preparation
# Load data and scale continuous covariates
dt = read_spr_data("../data/adbc23_spr.csv");

# exclude data
before = nrow(combine(groupby(dt, [:Condition, :Item, :Subject]), [:ReadingTime => mean => :ReadingTime]));
dt = exclude_trial(dt[((dt.Region .!= "Pre-critical_2")),:], 50, 2500, 50, 6000);
after = nrow(combine(groupby(dt, [:Condition, :Item, :Subject]), [:ReadingTime => mean => :ReadingTime]))
round((before - after) / before * 100, digits=2)

## Model Preparation
contr = Dict(:Subject => Grouping(), :Item => Grouping())
f1 = @formula(logRT ~ 1 + Plaus + Cloze_distractor + (1 + Plaus + Cloze_distractor | Item) + (1 + Plaus + Cloze_distractor | Subject))
f2 = @formula(logRT ~ 1 + PrecritRT + Plaus + Cloze_distractor + (1 + PrecritRT + Plaus + Cloze_distractor | Item) + (1 +  PrecritRT + Plaus + Cloze_distractor | Subject))

## Run Models
f = f1
lmerSPR=fit_models(dt, f, contr);
numpred = length(f.rhs) - 2
coefs = [Symbol(x) for x in string.(f.rhs[1:numpred])];
lmerSPR = combine_datasets(dt, lmerSPR, coefs[2:length(coefs)])
lmerSPR = compute_RTs(lmerSPR, coefs);
ests = names(lmerSPR)[occursin.("est_", names(lmerSPR))]
lmerSPR = compute_residuals(lmerSPR, ests); 
lmerSPR = compute_coefs(lmerSPR, coefs); 

CSV.write("../data/rRTs_Plaus_Clozedist.csv", lmerSPR);

f = f2
lmerSPR = fit_models(dt, f, contr);
numpred = length(f.rhs) - 2
coefs = [Symbol(x) for x in string.(f.rhs[1:numpred])];
lmerSPR = combine_datasets(dt, lmerSPR, coefs[2:length(coefs)])
lmerSPR = compute_RTs(lmerSPR, coefs);
ests = names(lmerSPR)[occursin.("est_", names(lmerSPR))]
lmerSPR = compute_residuals(lmerSPR, ests); 
lmerSPR = compute_coefs(lmerSPR, coefs); 

CSV.write("../data/rRTs_PrecritRT_Plaus_Clozedist.csv", lmerSPR);