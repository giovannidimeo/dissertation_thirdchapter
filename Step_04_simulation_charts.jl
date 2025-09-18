# Can be run only after running step 2 and 3 
# Remember to load directories and data from step 2

using Measures
######################################################################

df_main = copy(dfs[1])
df_counterfactual1 = copy(dfs[2])

######################################################################
# APPEND DataFrames
df_main.counterfactual .= 0 
df_counterfactual1.counterfactual .= 1
df_counterfactual1.id .= df_counterfactual1.id .+ 10000

df = vcat(df_main,df_counterfactual1) 

df.h_costs = df.HT .* df.pt .* (ψ_r +  δ + sell_costs1) .* df.DP .* df.DT  + 
df.HT .* df.pt .* (1 .- df.OT) .* (ψ_r) .* (1 .- df.DP) .* df.DT .* (1 .- df.DS) +
df.HT .* df.pt .* (1 .- df.OT) .* (ψ_r + sell_costs) .* df.DT .* (df.DS) + 
df.HT .* df.pt .* (ψ_r) .* (1 .- df.DT)

df.tot_wealth = df.ST .+ df.DT .* df.HT .* df.OT .* df.pt

df.event .= df.time .- df.SH

######################################################################
# CONVERT PRICES 
cd(background_data)
index_df = CSV.read("index.csv", DataFrame)
index_df.time .= index_df.year .- 1977
select!(index_df, [:time, :year, :index_forecast])
leftjoin!(df, index_df, on = :time) 

df.ST .= df.ST ./ df.index_forecast
df.tot_wealth .= df.tot_wealth ./ df.index_forecast
df.h_costs .= df.h_costs ./ df.index_forecast
df.CT .= df.CT ./ df.index_forecast
df.YT .= df.YT ./ df.index_forecast

df.tot_consumption .= (df.CT .+ df.h_costs)

df.ST = winsor2(df.ST/10, 95)
df.tot_wealth = winsor2(df.tot_wealth/10, 95)
df.h_costs =  df.h_costs/10
df.CT =  winsor2(df.CT/10, 95)
df.tot_consumption .= winsor2(df.tot_consumption/10, 95)


cd(results)
CSV.write("dataready.csv", df)


default(;fontfamily="serif-roman")

######################################################################
# EVENT STUDIES 
function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event, :counterfactual]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end


# Prob of purchasing
df_mean = collapsed(df, "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "B: House purchase",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p2 = plot!()

mean(@rsubset(df_mean, :counterfactual == 0, :event >= 0).avg)

mean(@rsubset(df_mean, :counterfactual == 1, :event >= 0).avg)


# Ownership rate
df_mean = collapsed(df, "DT")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "A: Homeownership",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)


plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)


plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p1 = plot!()


@rsubset(df_mean, :counterfactual == 1)
@rsubset(df_mean, :counterfactual == 0)


######################################################################
# Survival analysis 

cd(results)
survival_data = CSV.read("survival.csv", DataFrame)

plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Treatment",
title = "C: Survival function",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data, :group == 0).time, 
@rsubset(survival_data,  :group == 0).surv,
ribbon = (@rsubset(survival_data, :group == 0).upper .- @rsubset(survival_data, :group == 0).surv, 
@rsubset(survival_data, :group == 0).surv .- @rsubset(survival_data, :group == 0).lower),
label="Baseline",
linestyle = :dashdot)

p3 = plot!()

@rsubset(survival_data, :group == 0, :time == 30)
@rsubset(survival_data, :group == 1, :time == 30)

######################################################################
# Non-housing consumption

df.owner .= ifelse.(df.event .== 0 .&& df.DT .== 1,  1, 0)

transform!(groupby(df, [:id]), :owner => maximum => :owner)

df_restricted = select(df, [:id, :event, :CT, :counterfactual, :owner])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted = unstack(df_restricted, [:id, :event], :counterfactual, :CT, renamecols = x -> Symbol("Group", x))

leftjoin!(df_restricted, select(df, [:id, :event, :owner]), on = [:id, :event]) 

df_restricted.ratio .= df_restricted.Group1 ./ df_restricted.Group0
df_restricted.ratio .= ifelse.(isnan.(df_restricted.ratio) ,  1, df_restricted.ratio)

df_owners = @rsubset(df_restricted,  :owner == 1)
df_tenants = @rsubset(df_restricted,  :owner == 0)


function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

df_mean_t = collapsed(df_tenants, "ratio")
df_mean_o = collapsed(df_owners, "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "D: Non-housing consumption",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants",
linestyle = :dashdot)

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
savefig("purchases_sales.pdf")


######################################################################
# Age differences

function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event, :counterfactual]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

mean(@rsubset(df, :SH <= 20, :time == 40, :counterfactual == 1).DT)
mean(@rsubset(df, :SH > 20, :time == 40, :counterfactual == 1).DT)

mean(@rsubset(df, :SH <= 20, :time == 40, :counterfactual == 0).DT)
mean(@rsubset(df, :SH > 20, :time == 40, :counterfactual == 0).DT)


# Prob of purchasing
df_mean = collapsed(@rsubset(df, :SH <= 20), "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "A: House purchase - younger",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p1 = plot!()


df_mean = collapsed(@rsubset(df, :SH > 20), "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "B: House purchase - older",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p2 = plot!()


# Survival analysis
cd(results)
survival_data = CSV.read("survival_younger.csv", DataFrame)

plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Treatment",
title = "C: Survival function - younger",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data, :group == 0).time, 
@rsubset(survival_data,  :group == 0).surv,
ribbon = (@rsubset(survival_data, :group == 0).upper .- @rsubset(survival_data, :group == 0).surv, 
@rsubset(survival_data, :group == 0).surv .- @rsubset(survival_data, :group == 0).lower),
label="Baseline",
linestyle = :dashdot)

p3 = plot!()

survival_data = CSV.read("survival_older.csv", DataFrame)

plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Treatment",
title = "D: Survival function - older",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data, :group == 0).time, 
@rsubset(survival_data,  :group == 0).surv,
ribbon = (@rsubset(survival_data, :group == 0).upper .- @rsubset(survival_data, :group == 0).surv, 
@rsubset(survival_data, :group == 0).surv .- @rsubset(survival_data, :group == 0).lower),
label="Baseline",
linestyle = :dashdot)

p4 = plot!()

# Non-housing consumption
df.owner .= ifelse.(df.event .== 0 .&& df.DT .== 1,  1, 0)

transform!(groupby(df, [:id]), :owner => maximum => :owner)

df_restricted = select(df, [:id, :event, :CT, :counterfactual, :owner])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted = unstack(df_restricted, [:id, :event], :counterfactual, :CT, renamecols = x -> Symbol("Group", x))

leftjoin!(df_restricted, select(df, [:id, :event, :owner, :SH]), on = [:id, :event]) 

df_restricted.ratio .= df_restricted.Group1 ./ df_restricted.Group0
df_restricted.ratio .= ifelse.(isnan.(df_restricted.ratio) ,  1, df_restricted.ratio)

df_owners = @rsubset(df_restricted,  :owner == 1)
df_tenants = @rsubset(df_restricted,  :owner == 0)


function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

df_mean_t = collapsed(@rsubset(df_tenants, :SH <= 20), "ratio")
df_mean_o = collapsed(@rsubset(df_owners, :SH <= 20), "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "E: Non-housing consumption - younger",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p5 = plot!()

df_mean_t = collapsed(@rsubset(df_tenants, :SH > 20), "ratio")
df_mean_o = collapsed(@rsubset(df_owners, :SH > 20), "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "F: Non-housing consumption - older",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p6 = plot!()

cd(results)
plot(p1, p2, p3, p4, p5, p6, layout = grid(3,2))
plot!(size=(1000,1500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig("age_het.pdf")


######################################################################
# Income differences
df.income_flag .= ifelse.(df.event .== -1 .&& df.YT .<= median(@rsubset(df, :event == -1).YT), 0, 1)
df.income_flag .= ifelse.(df.event .!= -1, 99, df.income_flag)
transform!(groupby(df, [:id]), :income_flag => minimum => :income_flag)


function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event, :counterfactual]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

mean(@rsubset(df, :income_flag == 0, :time == 40, :counterfactual == 0).DT)
mean(@rsubset(df, :income_flag== 0, :time == 40, :counterfactual == 1).DT)

mean(@rsubset(df, :income_flag == 1, :time == 40, :counterfactual == 0).DT)
mean(@rsubset(df, :income_flag == 1, :time == 40, :counterfactual == 1).DT)


# Prob of purchasing
df_mean = collapsed(@rsubset(df, :income_flag == 0), "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "A: House purchase - low income",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p1 = plot!()

df_mean = collapsed(@rsubset(df, :income_flag == 1), "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "B: House purchase - high income",
xlabel = "Years from shock",
ylabel = "Share",
label="Treatment",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean, :counterfactual == 0).event, 
@rsubset(df_mean,  :counterfactual == 0).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 0).hi .- @rsubset(df_mean, :counterfactual == 0).avg, 
@rsubset(df_mean, :counterfactual == 0).avg .- @rsubset(df_mean, :counterfactual == 0).lo),
label="Baseline",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p2 = plot!()

# Prob of selling
survival_data = CSV.read("survival_poorer.csv", DataFrame)

@rsubset(survival_data, :group == 0)
@rsubset(survival_data, :group == 1)


plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Treatment",
title = "C: Survival function - low income",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data, :group == 0).time, 
@rsubset(survival_data,  :group == 0).surv,
ribbon = (@rsubset(survival_data, :group == 0).upper .- @rsubset(survival_data, :group == 0).surv, 
@rsubset(survival_data, :group == 0).surv .- @rsubset(survival_data, :group == 0).lower),
label="Baseline",
linestyle = :dashdot)

p3 = plot!()

survival_data = CSV.read("survival_richer.csv", DataFrame)

plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Treatment",
title = "D: Survival function - high income",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data, :group == 0).time, 
@rsubset(survival_data,  :group == 0).surv,
ribbon = (@rsubset(survival_data, :group == 0).upper .- @rsubset(survival_data, :group == 0).surv, 
@rsubset(survival_data, :group == 0).surv .- @rsubset(survival_data, :group == 0).lower),
label="Baseline",
linestyle = :dashdot)

p4 = plot!()

# Non-housing consumption
df.owner .= ifelse.(df.event .== 0 .&& df.DT .== 1,  1, 0)
transform!(groupby(df, [:id]), :owner => maximum => :owner)
df_restricted = select(df, [:id, :event, :CT, :counterfactual, :owner])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted = unstack(df_restricted, [:id, :event], :counterfactual, :CT, renamecols = x -> Symbol("Group", x))

leftjoin!(df_restricted, select(df, [:id, :event, :owner, :income_flag]), on = [:id, :event]) 

df_restricted.ratio .= df_restricted.Group1 ./ df_restricted.Group0
df_restricted.ratio .= ifelse.(isnan.(df_restricted.ratio) ,  1, df_restricted.ratio)

df_owners = @rsubset(df_restricted,  :owner == 1)
df_tenants = @rsubset(df_restricted,  :owner == 0)


function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

df_mean_t = collapsed(@rsubset(df_tenants, :income_flag == 0), "ratio")
df_mean_o = collapsed(@rsubset(df_owners, :income_flag == 0), "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "E: Non-housing consumption - low income",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p5 = plot!()

df_mean_t = collapsed(@rsubset(df_tenants, :income_flag == 1), "ratio")
df_mean_o = collapsed(@rsubset(df_owners, :income_flag == 1), "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "F: Non-housing consumption - high income",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p6 = plot!()

cd(results)
plot(p1, p2, p3, p4, p5, p6, layout = grid(3,2))
plot!(size=(1000,1500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig("income_het.pdf")


### 
# Analysis of wealth
df.owner .= ifelse.(df.event .== 0 .&& df.DT .== 1,  1, 0)
transform!(groupby(df, [:id]), :owner => maximum => :owner)
df_restricted = select(df, [:id, :event, :tot_wealth, :counterfactual, :owner])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted = unstack(df_restricted, [:id, :event], :counterfactual, :tot_wealth, renamecols = x -> Symbol("Group", x))
leftjoin!(df_restricted, select(df, [:id, :event, :time, :owner, :income_flag]), on = [:id, :event]) 
df_restricted.ratio .= df_restricted.Group1 ./ df_restricted.Group0
df_restricted.ratio .= ifelse.(isnan.(df_restricted.ratio) ,  1, df_restricted.ratio)
df_restricted = @rsubset(df_restricted,  :ratio < 10) # get rid of outliers

df_owners = @rsubset(df_restricted,  :owner == 1)
df_restricted = @rsubset(df_restricted,  :owner == 0)

function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end

df_mean_t = collapsed(@rsubset(df_tenants, :income_flag == 1), "ratio")
df_mean_o = collapsed(@rsubset(df_owners, :income_flag == 1), "ratio")

plot(@rsubset(df_mean_o).event, 
@rsubset(df_mean_o).avg,
ribbon = (@rsubset(df_mean_o).hi .- @rsubset(df_mean_o).avg, 
@rsubset(df_mean_o).avg .- @rsubset(df_mean_o).lo),
title = "Total net wealth - high income",
xlabel = "Years from shock",
ylabel = "Ratio",
label="Owners",
xticks = -5:5:30,
ylimits = (0.5,1),
yticks = 0.5:0.1:1,
color = "green",
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(df_mean_t).event, 
@rsubset(df_mean_t).avg,
ribbon = (@rsubset(df_mean_t).hi .- @rsubset(df_mean_t).avg, 
@rsubset(df_mean_t).avg .- @rsubset(df_mean_t).lo),
color = "#F5E81F",
label="Tenants")

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "black", 
label = "")


######################################################################
# Elasticity of homeownership
df.income_flag .= ifelse.(df.event .== -1 .&& df.YT .<= median(@rsubset(df, :event == -1).YT), 0, 1)
df.income_flag .= ifelse.(df.event .!= -1, 99, df.income_flag)
transform!(groupby(df, [:id]), :income_flag => minimum => :income_flag)

df_restricted = select(df, [:id, :event, :DT, :YT, :HT, :SH, :income_flag, :counterfactual])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted.HT .= ifelse.(df_restricted.DT .== 0 , 0, df_restricted.HT)

# - All
df_restricted1 = unstack(df_restricted, [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(df_restricted, [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean = @rsubset(df_mean, :event >= 0)
df_mean = @rsubset(df_mean, :event <= 30)
df_mean.change .= ((df_mean.DT1_mean .- df_mean.DT0_mean) ./ (df_mean.DT0_mean)) ./ ((df_mean.YT1_mean .- df_mean.YT0_mean) ./ (df_mean.YT0_mean))

mean(df_mean.change)

# - Young
df_restricted1 = unstack(@rsubset(df_restricted, :SH <= 20), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :SH <= 20), [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_y = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean_y = @rsubset(df_mean_y, :event >= 0)
df_mean_y = @rsubset(df_mean_y, :event <= 30)
df_mean_y.change .= ((df_mean_y.DT1_mean .- df_mean_y.DT0_mean) ./ (df_mean_y.DT0_mean)) ./ ((df_mean_y.YT1_mean .- df_mean_y.YT0_mean) ./ (df_mean_y.YT0_mean))

# - Old
df_restricted1 = unstack(@rsubset(df_restricted, :SH > 20), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :SH > 20), [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_o = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean_o = @rsubset(df_mean_o, :event >= 0)
df_mean_o = @rsubset(df_mean_o, :event <= 30)
df_mean_o.change .= ((df_mean_o.DT1_mean .- df_mean_o.DT0_mean) ./ (df_mean_o.DT0_mean)) ./ ((df_mean_o.YT1_mean .- df_mean_o.YT0_mean) ./ (df_mean_o.YT0_mean))

# - Poor
df_restricted1 = unstack(@rsubset(df_restricted, :income_flag == 0), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :income_flag == 0), [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_p = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean_p = @rsubset(df_mean_p, :event >= 0)
df_mean_p = @rsubset(df_mean_p, :event <= 30)
df_mean_p.change .= ((df_mean_p.DT1_mean .- df_mean_p.DT0_mean) ./ (df_mean_p.DT0_mean)) ./ ((df_mean_p.YT1_mean .- df_mean_p.YT0_mean) ./ (df_mean_p.YT0_mean))

# - Rich 
df_restricted1 = unstack(@rsubset(df_restricted, :income_flag == 1), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :income_flag == 1), [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_r = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean_r = @rsubset(df_mean_r, :event >= 0)
df_mean_r = @rsubset(df_mean_r, :event <= 30)
df_mean_r.change .= ((df_mean_r.DT1_mean .- df_mean_r.DT0_mean) ./ (df_mean_r.DT0_mean)) ./ ((df_mean_r.YT1_mean .- df_mean_r.YT0_mean) ./ (df_mean_r.YT0_mean))


plot(df_mean.event, 
df_mean.change,
title = "A: Elasticity of homeownership",
xlabel = "Years from shock",
ylabel = "Elasticity",
xticks = 0:5:30,
ylimits = (0,2.5),
yticks = 0:0.5:2.5,
label = "All",
legend = :topleft,
legend_font_pointsize = 11)

plot!(df_mean_y.event, 
df_mean_y.change,
label = "Younger",
linestyle = :dash)

plot!(df_mean_o.event, 
df_mean_o.change,
label = "Older",
linestyle = :dot)

plot!(df_mean_p.event, 
df_mean_p.change,
label = "Low income",
linestyle = :dashdot)

plot!(df_mean_r.event, 
df_mean_r.change,
label = "High income",
linestyle = :dashdotdot)

p1 = plot!()


# Elasticity of housing size
# - All
df_restricted1 = unstack(df_restricted, [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(df_restricted, [:id, :event], :counterfactual, :HT, renamecols = x -> Symbol("HT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :HT0, :HT1] .=> mean) 
df_mean = @rsubset(df_mean, :event >= 0)
df_mean = @rsubset(df_mean, :event <= 30)
df_mean.change .= ((df_mean.HT1_mean .- df_mean.HT0_mean) ./ (df_mean.HT0_mean)) ./ ((df_mean.YT1_mean .- df_mean.YT0_mean) ./ (df_mean.YT0_mean))

mean(df_mean.change)

# - Young
df_restricted1 = unstack(@rsubset(df_restricted, :SH <= 20), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :SH <= 20), [:id, :event], :counterfactual, :HT, renamecols = x -> Symbol("HT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_y = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :HT0, :HT1] .=> mean) 
df_mean_y = @rsubset(df_mean_y, :event >= 0)
df_mean_y = @rsubset(df_mean_y, :event <= 30)
df_mean_y.change .= ((df_mean_y.HT1_mean .- df_mean_y.HT0_mean) ./ (df_mean_y.HT0_mean)) ./ ((df_mean_y.YT1_mean .- df_mean_y.YT0_mean) ./ (df_mean_y.YT0_mean))

mean(df_mean.change)

# - Old
df_restricted1 = unstack(@rsubset(df_restricted, :SH > 20), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :SH > 20), [:id, :event], :counterfactual, :HT, renamecols = x -> Symbol("HT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_o = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :HT0, :HT1] .=> mean) 
df_mean_o = @rsubset(df_mean_o, :event >= 0)
df_mean_o = @rsubset(df_mean_o, :event <= 30)
df_mean_o.change .= ((df_mean_o.HT1_mean .- df_mean_o.HT0_mean) ./ (df_mean_o.HT0_mean)) ./ ((df_mean_o.YT1_mean .- df_mean_o.YT0_mean) ./ (df_mean_o.YT0_mean))

# - Poor
df_restricted1 = unstack(@rsubset(df_restricted, :income_flag == 0), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :income_flag == 0), [:id, :event], :counterfactual, :HT, renamecols = x -> Symbol("HT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_p = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :HT0, :HT1] .=> mean) 
df_mean_p = @rsubset(df_mean_p, :event >= 0)
df_mean_p = @rsubset(df_mean_p, :event <= 30)
df_mean_p.change .= ((df_mean_p.HT1_mean .- df_mean_p.HT0_mean) ./ (df_mean_p.HT0_mean)) ./ ((df_mean_p.YT1_mean .- df_mean_p.YT0_mean) ./ (df_mean_p.YT0_mean))

# - Rich 
df_restricted1 = unstack(@rsubset(df_restricted, :income_flag == 1), [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(@rsubset(df_restricted, :income_flag == 1), [:id, :event], :counterfactual, :HT, renamecols = x -> Symbol("HT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean_r = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :HT0, :HT1] .=> mean) 
df_mean_r = @rsubset(df_mean_r, :event >= 0)
df_mean_r = @rsubset(df_mean_r, :event <= 30)
df_mean_r.change .= ((df_mean_r.HT1_mean .- df_mean_r.HT0_mean) ./ (df_mean_r.HT0_mean)) ./ ((df_mean_r.YT1_mean .- df_mean_r.YT0_mean) ./ (df_mean_r.YT0_mean))


plot(df_mean.event, 
df_mean.change,
title = "B: Elasticity of owned housing size",
xlabel = "Years from shock",
ylabel = "Elasticity",
xticks = 0:5:30,
ylimits = (0,2.5),
yticks = 0:0.5:2.5,
label = "All",
legend = :topleft,
legend_font_pointsize = 11)

plot!(df_mean_y.event, 
df_mean_y.change,
label = "Younger",
linestyle = :dash)

plot!(df_mean_o.event, 
df_mean_o.change,
label = "Older",
linestyle = :dot)

plot!(df_mean_p.event, 
df_mean_p.change,
label = "Low income",
linestyle = :dashdot)

plot!(df_mean_r.event, 
df_mean_r.change,
label = "High income",
linestyle = :dashdotdot)

p2 = plot!()

plot(p1, p2, layour = grid(2,1))
plot!(size=(1000,500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
cd(results)
savefig("elasticity.pdf")




######################################################################
#  No transaction costs

df_main_app.counterfactual .= 0 
df_counterfactual1_app.counterfactual .= 1
df_counterfactual1_app.id .= df_counterfactual1.id .+ 10000

df_app = vcat(df_main_app,df_counterfactual1_app) 
df_app.event .= df_app.time .- df_app.SH

cd(background_data)
index_df = CSV.read("index.csv", DataFrame)
index_df.time .= index_df.year .- 1977
select!(index_df, [:time, :year, :index_forecast])
leftjoin!(df_app, index_df, on = :time) 

df_app.YT .= df_app.YT ./ df_app.index_forecast

cd(results)
CSV.write("dataready_app.csv", df_app)


function collapsed(df, var)
    var_sym = Symbol(var)
    df_mean = combine(groupby(df, [:event, :counterfactual]), var_sym .=> mean .=> "avg", var_sym .=> std .=> :se, nrow => :count, renamecols = false)
    df_mean = @rsubset(df_mean, :event >= -5, :event <= 30)
    
    df_mean.hi = df_mean.avg .+ 1.96 .* df_mean.se ./ sqrt.(df_mean.count)
    df_mean.lo = df_mean.avg .- 1.96 .* df_mean.se ./ sqrt.(df_mean.count)

    return(df_mean)
end


# Homeownership
df_mean = collapsed(df, "DT")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "A: Homeownership",
xlabel = "Years from shock",
ylabel = "Share",
label="Transaction costs",
xticks = -5:5:30,
ylimits = (0,1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

df_mean_app = collapsed(df_app, "DT")

plot!(@rsubset(df_mean_app, :counterfactual == 1).event, 
@rsubset(df_mean_app, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean_app, :counterfactual == 1).hi .- @rsubset(df_mean_app, :counterfactual == 1).avg, 
@rsubset(df_mean_app, :counterfactual == 1).avg .- @rsubset(df_mean_app, :counterfactual == 1).lo),
label="No transaction costs",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p1 = plot!()

# Purchases
df_mean = collapsed(df, "DP")

plot(@rsubset(df_mean, :counterfactual == 1).event, 
@rsubset(df_mean, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean, :counterfactual == 1).hi .- @rsubset(df_mean, :counterfactual == 1).avg, 
@rsubset(df_mean, :counterfactual == 1).avg .- @rsubset(df_mean, :counterfactual == 1).lo),
title = "B: House purchase",
xlabel = "Years from shock",
ylabel = "Share",
label="Transaction costs",
xticks = -5:5:30,
ylimits = (0,0.1),
yticks = 0:0.02:0.1,
legend = :topright,
legend_font_pointsize = 11)

df_mean_app = collapsed(df_app, "DP")

plot!(@rsubset(df_mean_app, :counterfactual == 1).event, 
@rsubset(df_mean_app, :counterfactual == 1).avg,
ribbon = (@rsubset(df_mean_app, :counterfactual == 1).hi .- @rsubset(df_mean_app, :counterfactual == 1).avg, 
@rsubset(df_mean_app, :counterfactual == 1).avg .- @rsubset(df_mean_app, :counterfactual == 1).lo),
label="No transaction costs",
linestyle = :dashdot)

plot!([0], 
seriestype = "vline", 
linestyle = :dash, 
linecolor = "red", 
label = "")

p2 = plot!()


# Elasticity of homeownership
df_restricted = select(df, [:id, :event, :DT, :YT, :HT, :ST, :counterfactual])
df_restricted.id .= df_restricted.id .- (df_restricted.counterfactual .* 10000)
df_restricted.HT .= ifelse.(df_restricted.DT .== 0 , 0, df_restricted.HT)

df_restricted1 = unstack(df_restricted, [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_restricted2 = unstack(df_restricted, [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_restricted1 = innerjoin(df_restricted1, df_restricted2, on = [:id, :event])
df_mean = combine(groupby(df_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_mean = @rsubset(df_mean, :event >= 0)
df_mean = @rsubset(df_mean, :event <= 30)
df_mean.change .= ((df_mean.DT1_mean .- df_mean.DT0_mean) ./ (df_mean.DT0_mean)) ./ ((df_mean.YT1_mean .- df_mean.YT0_mean) ./ (df_mean.YT0_mean))

df_app_restricted = select(df_app, [:id, :event, :DT, :YT, :HT, :SH, :counterfactual])
df_app_restricted.id .= df_app_restricted.id .- (df_app_restricted.counterfactual .* 10000)
df_app_restricted.HT .= ifelse.(df_app_restricted.DT .== 0 , 0, df_app_restricted.HT)

df_app_restricted1 = unstack(df_app_restricted, [:id, :event], :counterfactual, :YT, renamecols = x -> Symbol("YT", x))
df_app_restricted2 = unstack(df_app_restricted, [:id, :event], :counterfactual, :DT, renamecols = x -> Symbol("DT", x))
df_app_restricted1 = innerjoin(df_app_restricted1, df_app_restricted2, on = [:id, :event])
df_app_mean = combine(groupby(df_app_restricted1, [:event]), [:YT0, :YT1, :DT0, :DT1] .=> mean) 
df_app_mean = @rsubset(df_app_mean, :event >= 0)
df_app_mean = @rsubset(df_app_mean, :event <= 30)
df_app_mean.change .= ((df_app_mean.DT1_mean .- df_app_mean.DT0_mean) ./ (df_app_mean.DT0_mean)) ./ ((df_app_mean.YT1_mean .- df_app_mean.YT0_mean) ./ (df_app_mean.YT0_mean))


plot(df_mean.event, 
df_mean.change,
title = "D: Elasticity of homeownership",
xlabel = "Years from shock",
ylabel = "Elasticity",
xticks = 0:5:30,
ylimits = (0,2.5),
yticks = 0:0.5:2.5,
label = "Transaction costs",
legend = :topright,
legend_font_pointsize = 11)

plot!(df_app_mean.event, 
df_app_mean.change,
label = "No transaction costs",
linestyle = :dashdot)

p4 = plot!()


cd(results)
survival_data = CSV.read("survival.csv", DataFrame)
survival_data_app = CSV.read("survival_App.csv", DataFrame)

plot(@rsubset(survival_data, :group == 1).time, 
@rsubset(survival_data,  :group == 1).surv,
ribbon = (@rsubset(survival_data, :group == 1).upper .- @rsubset(survival_data, :group == 1).surv, 
@rsubset(survival_data, :group == 1).surv .- @rsubset(survival_data, :group == 1).lower),
label="Transaction costs",
title = "C: Survival function",
xlabel = "Years from shock",
ylabel = "Survival rate",
xticks = 0:5:30,
ylimits = (0, 1),
yticks = 0:0.2:1,
legend = :topright,
legend_font_pointsize = 11)

plot!(@rsubset(survival_data_app, :group == 1).time, 
@rsubset(survival_data_app,  :group == 1).surv,
ribbon = (@rsubset(survival_data_app, :group == 1).upper .- @rsubset(survival_data_app, :group == 1).surv, 
@rsubset(survival_data_app, :group == 1).surv .- @rsubset(survival_data_app, :group == 1).lower),
label="No transaction costs",
linestyle = :dashdot)

p3 = plot!()



plot(p1, p2, p3, p4, layour = grid(2,2))
plot!(size=(1000,1000), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
cd(results)
savefig("transaction_costs.pdf")
