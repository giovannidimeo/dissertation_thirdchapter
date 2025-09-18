# Can be run only after running step 2 and 5
# Remember to load directories and data from step 2

using StatsPlots
using Measures

######################################################################
# SHP, descriptive stats 

# Histogram
cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))
shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
shp = @rsubset(shp, :DT > -1)
shp = @rsubset(shp, ((:shock .> 0) .& (:increase_in_neg1 .== 0)))
shp.age_shock  .= shp.age - shp.event
shp = @rsubset(shp, :max_age <= 65)
shp = @rsubset(shp, :age_shock <= 65)

shp.data  .= 1
shp.event .= shp.time .- shp.shock
transform!(groupby(shp, [:idhous]), :event => minimum  => :event_min)
transform!(groupby(shp, [:idhous]), eachindex => :idnr)
shp = @rsubset(shp, :idnr == 1) #take only one observartion per id

df_cat = select(shp, [:id, :time, :age, :DT, :YT, :data, :shock])

df_cat.age_cat .= 0
df_cat.age_cat.= ifelse.((df_cat.age .<= 25), 1, df_cat.age_cat )
for l = 1:20
    q1 = 25 + 2*(l-1)
    q2 = 25 + 2*l
    # Use df_cat.age instead of df_cat.age_cat in the condition
    df_cat.age_cat .= ifelse.((df_cat.age .> q1) .& (df_cat.age .<= q2), l+1, df_cat.age_cat)
end

df_mean_tot = combine(groupby(df_cat, [:data]), nrow => :counts, renamecols = true)
df_mean = combine(groupby(df_cat, [:data, :age_cat]), nrow => :tot, renamecols = false)
leftjoin!(df_mean, df_mean_tot, on = :data) 
df_mean.ratio .= df_mean.tot ./ df_mean.counts

df_mean.age_cat2 .= df_mean.age_cat .+ (df_mean.age_cat .-1) .+ 24

default(;fontfamily="serif-roman")
cd(background_data)
p1 = groupedbar(df_mean.age_cat2, df_mean.ratio, group = df_mean.data, 
yticks = (0:0.02:0.1),
ylimits = (0,0.1),
grid = :on,
lw = 0,
xticks = (25:5:65),
title = "A: Average Household age at shock",
xlabel = "Age",
ylabel = "Frequency",
legend = false
)

# Income 
cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))
shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
shp = @rsubset(shp, :DT > -1)
shp = @rsubset(shp, ((:shock .> 0) .& (:increase_in_neg1 .== 0)))
shp.age_shock  .= shp.age - shp.event
shp = @rsubset(shp, :max_age <= 65)
shp = @rsubset(shp, :age_shock <= 65)
shp.event .= shp.time .- shp.shock
transform!(groupby(shp, [:idhous]), :event => minimum  => :event_min)

shp_shock = copy(shp)

select!(shp_shock, [:id, :time, :age, :DT, :YT])

cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))
shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
transform!(groupby(shp, [:idhous]), :max_age => minimum  => :event_min)
shp = @rsubset(shp, :max_age < 60)
shp = @rsubset(shp, :DT > -1)
shp = @rsubset(shp, (:shock .== 0))

shp_control= copy(shp)

select!(shp_control, [:id, :time, :age, :DT, :YT])

lm1 = lm(@formula(YT ~ age + age^2  +  age^3 ), shp_shock)
ages = collect(25:1:65)
newdata_shp = DataFrame(age = ages)
newdata_shp.age2 = newdata_shp.age .^ 2
newdata_shp.age3 = newdata_shp.age .^ 3

predicted_YT_shp = predict(lm1, newdata_shp, interval=:confidence)

# Polynomial in simulated data
lm1 = lm(@formula(YT ~ age + age^2  +  age^3 ), shp_control)
ages = collect(25:1:65)
newdata_simul = DataFrame(age = ages)
newdata_simul.age2 = newdata_simul.age .^ 2
newdata_simul.age3 = newdata_simul.age .^ 3

predicted_YT_simul = predict(lm1, newdata_simul, interval=:confidence)

plot(ages, 
predicted_YT_simul.prediction,
ribbon = (predicted_YT_simul.upper - predicted_YT_simul.prediction, predicted_YT_simul.prediction - predicted_YT_simul.lower),
title = "B: Per-capita household income",
xlabel = "Age",
ylabel = "100,000 CHF",
label="Treated",
xticks = 25:5:65,
ylimits = (0,1),
legend = :topleft,
legend_font_pointsize = 11)

plot!(ages, 
predicted_YT_shp.prediction,
ribbon = (predicted_YT_shp.upper - predicted_YT_shp.prediction, predicted_YT_shp.prediction - predicted_YT_shp.lower),
label="Control",
linestyle = :dashdot)

p2 =plot!()

# Evolution of income by event time
cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))
shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
shp = @rsubset(shp, :DT > -1)
shp = @rsubset(shp, ((:shock .> 0) .& (:increase_in_neg1 .== 0)))
shp.age_shock  .= shp.age - shp.event
shp = @rsubset(shp, :max_age <= 65)
shp = @rsubset(shp, :age_shock <= 65)
shp.event .= shp.time .- shp.shock

function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 5)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

#INCOME 
df_mean = collapsed(shp, "YT")

sort!(df_mean, [:event])

plot(df_mean.event, 
df_mean.avg,
ribbon = (df_mean.hi .- df_mean.avg, 
df_mean.avg .- df_mean.lo),
title = "C: Evolution of p.c. household income",
xlabel = "Years from shock",
ylabel = "100,000 CHF",
#linecolor = "black",
xticks = -5:1:10,
ylimits = (0,1),
legend = false)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p3 = plot!()


#HOMEOWNERSHIP 
df_mean = collapsed(shp, "DT")
sort!(df_mean, [:event])

plot(df_mean.event, 
df_mean.avg,
ribbon = (df_mean.hi .- df_mean.avg, 
df_mean.avg .- df_mean.lo),
title = "D: Evolution of homeownership",
xlabel = "Years from shock",
ylabel = "Share",
#linecolor = "black",
xticks = -5:1:5,
ylimits = (0,1),
legend = false)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p4 = plot!()

cd(results)
plot(p1, p2, p3, p4, layour = (2,2))
plot!(size=(1000,1000), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig("desc_shp.pdf")

######################################################################
# SHP, empirical results 
cd(results)
rresults = load("AllResults_new.RData")
print("done")

for i = 1:4
    replace!(rresults["MainResults_Megaloop"][i].upper, missing => 0)
    replace!(rresults["MainResults_Megaloop"][i].lower, missing => 0)
end

plot(rresults["MainResults_Megaloop"]["results_simul_YT"].group, 
rresults["MainResults_Megaloop"]["results_simul_YT"].effect,
ribbon = (rresults["MainResults_Megaloop"]["results_simul_YT"].upper .- rresults["MainResults_Megaloop"]["results_simul_YT"].effect, 
rresults["MainResults_Megaloop"]["results_simul_YT"].effect .- rresults["MainResults_Megaloop"]["results_simul_YT"].lower),
linecolor = "#EA3680",
fillcolor = "#EA3680",
label="Simulation",
title = "A: Income",
xlabel = "Years from shock",
ylabel = "100,000 CHF",
xticks = -5:1:5,
ylimits = (-0.6,0.4),
yticks = -0.6:0.2:0.4,
legend = :topright,
legend_font_pointsize = 11)

plot!(rresults["MainResults_Megaloop"]["results_shp_YT"].group, 
rresults["MainResults_Megaloop"]["results_shp_YT"].effect,
ribbon = (rresults["MainResults_Megaloop"]["results_shp_YT"].upper .- rresults["MainResults_Megaloop"]["results_shp_YT"].effect, 
rresults["MainResults_Megaloop"]["results_shp_YT"].effect .- rresults["MainResults_Megaloop"]["results_shp_YT"].lower),
linecolor = "#000080",
fillcolor = "#000080",
label="SHP",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p1 = plot!()



plot(rresults["MainResults_Megaloop"]["results_simul_DT"].group, 
rresults["MainResults_Megaloop"]["results_simul_DT"].effect,
ribbon = (rresults["MainResults_Megaloop"]["results_simul_DT"].upper .- rresults["MainResults_Megaloop"]["results_simul_DT"].effect, 
rresults["MainResults_Megaloop"]["results_simul_DT"].effect .- rresults["MainResults_Megaloop"]["results_simul_DT"].lower),
linecolor = "#EA3680",
fillcolor = "#EA3680",
label="Simulation",
ylabel = "Share",
xticks = -5:1:5,
ylimits = (-0.3,0.2),
yticks = -0.3:0.1:0.2,
legend = :topright,
legend_font_pointsize = 11,
title = "B: Homeownership",
xlabel = "Years from shock" )

plot!(rresults["MainResults_Megaloop"]["results_shp_DT"].group, 
rresults["MainResults_Megaloop"]["results_shp_DT"].effect,
ribbon = (rresults["MainResults_Megaloop"]["results_shp_DT"].upper .- rresults["MainResults_Megaloop"]["results_shp_DT"].effect, 
rresults["MainResults_Megaloop"]["results_shp_DT"].effect .- rresults["MainResults_Megaloop"]["results_shp_DT"].lower),
linecolor = "#000080",
fillcolor = "#000080",
label="SHP",
linestyle = :dashdot) 

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p2= plot!()


cd(results)
plot(p1, p2, layour = (1,1))
plot!(size=(1000,500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig(raw"CS_shp_new.pdf")

####################
# DESCRIPTIVE ANALYSIS

df_counterfactual1 = dfs[2]
df_counterfactual1 = @rsubset(df_counterfactual1, :time < 40)
df_counterfactual1.data  .= 0
df_counterfactual1.shock .= df_counterfactual1.SH
df_counterfactual1.event .= df_counterfactual1.time .- df_counterfactual1.shock

function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 5)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

collapsed(df_counterfactual1, "DT")
