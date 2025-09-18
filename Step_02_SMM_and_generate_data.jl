begin
    using Pkg
    using BenchmarkTools
    using DataFrames
    using DataFramesMeta
    using StatsBase
    using Plots
    using XLSX
    using Distributions
    using Combinatorics
    using DataInterpolations
    using JLD
    using Roots
    using Random
    using Interpolations
    using CSV
    using Alert
    using GLM
    using Dates
end

####################################################################################################################
# Parameters and background data
begin # Set directories 
    background_data = raw"~\background_data"
    results = raw"~\results"
    results_smm = raw"~\calibration"
    shp_dir = raw"~" # folder where you stored SHP sample from stata
end 

begin
    points = 10
    T= 70  #length of life
    N= 40 #length of working life
    β = 0.97 # discont paramter
    γ = 5 # risk aversion parameter 
    δ = 0.2 # down payment
    η = 0.22 # preference for housing consumption 
    inv_eta = 1/(1-η)
    r_free = 0.018 # risk-free interest rate 
    R = 1 + r_free
    ϕ = 0.8 # replacement rate
    Random.seed!(42)
    sell_costs = 0.01 #0.01 # costs for selling 
    sell_costs1 = 0.05 #0.05 # costs for purchasing
end 

begin # Load data for housing prices
    cd(background_data)
    df_re = DataFrame(XLSX.readtable("re_indeces_prices.xlsx", "data"))

    df_re[!,:year] = string.(df_re[:,:year])
    df_re[!,:year] = parse.(Int, df_re[:,:year])

    ψ_r = mean(df_re.alpha_t) # Relative rental price

    # You can use the prices from the data. However, the series is very noisy. 
    pt = zeros(T) # Empty vector of housing prices
    Random.seed!(42)
    random = rand(Truncated(Normal(0,80), -10000, +10000), T)
    for t = 1:T
        pt[t] = 2463.248 + t*99.99899  + random[t] # Coefficients from simple OLS regression (using only years after 1977)
    end
end

begin # Create vector for mortality prob
    π1 = zeros(70)
    for x = 1:70
        π1[x] =  0.0002667*exp(0.1175*x)
    end

end

begin # Create income profile
    cd(background_data)
    y = DataFrame(XLSX.readtable("income.xlsx", "all_values"))

    y = y[:, [:year, :age_25, :age_35, :age_45, :age_55, :age_65]]
    y = @rsubset(y, :year > 1977)  # Drop 1977 because we don't have inflation data on this year
    y1 = Matrix(y[:, [:age_25, :age_35, :age_45, :age_55, :age_65]])
    rename!(y, [:"year", :"25", :"35", :"45", :"55", :"65"])
    y = stack(y, [:"25", :"35", :"45", :"55", :"65"], :year)
    rename!(y, [:variable, :value] .=> [:age, :avg_income])
    y[!, :year_standard] = y[!, :year] .- 1977
    y[!, :age_standard] = parse.(Float64, y[!, :age]) .- 24

    # Create interpolation for yearly income growth
    age_s = unique!(y[!, :age_standard])
    year_s = unique!(y[!, :year_standard])

    interp_y = extrapolate(interpolate((year_s,age_s,), y1, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation

    ypsilon = zeros(N)
    ypsilon[1] = 1

    for t in 2:N
        ypsilon[t] = 1 + (interp_y(t, t) - interp_y(t - 1, t - 1)) / interp_y(t - 1, t - 1) # yearly growth
    end

    # Unusual increase in the last years before retirment. Interpolate using only the first 35 years of life
    time_vect = range(1,N-5)
    interp_y_corrected = extrapolate(interpolate((time_vect,), ypsilon[1:N-5], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation

    for t in 2:N
        ypsilon[t] = interp_y_corrected(t)
    end

    mu = interp_y(1, 1) # mean income distribution for first year
    sd = interp_y(1, 1) * 0.6 # sd income distribution for first year
    lognormal_st = sqrt(log(1 + (sd/mu)^2)) # transform the values for lognormal distribution
    lognormal_mu = log(mu) - lognormal_st^2 / 2
 
    # Create possible income profiles
    first_income=  rand(LogNormal(lognormal_mu,lognormal_st), points) # draw  1000 starting incomes from lognormal distribution
    rndy= rand(Normal(0, 0.1,), (points,N)) 

    Y = zeros(Float64,(points, T)) # create matrix of income profiles
    Y[:,1] = first_income .* 12
    Y[:,1] = sort( Y[:,1])

    # Income profiles for working age
        for i in range(1, points)
            for t in range(1,N)
                Y[i, t+1] = (Y[i, t])*ypsilon[t]
            end
        end

    # Income profiles for retirement age
        for i in range(1, points)
            for t in range(N+1,T)
                Y[i,t] = Y[i,N] * ϕ
            end
        end
        
    # Sort values
    for t in range(1,T)
        Y[:,t] = sort( Y[:,t])
    end

end

begin # Create matrix for liquid wealth
    S = zeros(Float64,(points))
    S[1] = 0

    for i = 2:points
        S[i] = S[i-1]+(2500 + 8000*(i-2)^2)
        println(S[i])
    end
end

begin # Create matrix for net housing value
    ω = range(-0.1,0.8,points)
end 

begin # Create a vector with possible values of owned real esatate 
    h_bar = range(30,120, length = points)
end

begin #  Create sets of possible interest rates and income

    # Real estate
    rh_backw = [-1/100 2
    7/100 2
    -1/100 2
    7/100 2 ]

    # Income shocks
    downval = -0.3
    upval = 593/13600 
    prob = 4624/4673
    meanv = downval*(1-prob) + upval*prob
    ((downval-meanv)^2*(1-prob) + (upval-meanv)^2*prob)^0.5
    y_shocks_perm = [downval 1/(1-prob)
    upval 1/(prob)]

end

begin # Rescale vecotrs and matrices 
    S[:] .= S[:]/10000
    Y[:, :] .= Y[:, :]/10000
    pt[:] .= pt[:]/10000 
end

####################################################################################################################
# Define functions 

function h_c(ct, RE, τ) # Derive housing consumption for renters/owners
    if RE == 0
        ht = ct * (η/(1-η)) / (pt[T-τ] * ψ_r)
    else 
        ht = RE
    end
    return(ht)
end

function q(pval) # Total consumption costs for renters
    q = ((η/(1 - η))/(pval* ψ_r))^η
    return(q)
end

function u_c(ct, hval, pval) # Utility function
    if hval == 0 
        u_val = (ct * q(pval))^(1 - γ)/(1 - γ) 
    else 
        u_val = (ct^(1 - η) * hval^η)^(1 - γ)/(1 - γ)
    end
    return(u_val)
end

function v_w(w, ω_val, h, pval, k) # Utility of bequest  
    if h != 0
        v_w_val = (w + h_bar[h]*ω_val*pval + k)^(1-γ)/(1-γ)
    else
        v_w_val = (w + k)^(1-γ)/(1-γ)
    end
    return(v_w_val)
end

function second_upper_envelope_step(m_temp, c_temp, v_temp) # Second upper envelope
    j1 = Int[]
    j2 = Int[]
    j3 = Int[]
    vupper1 = Int[]
    vupper2 = Int[]
    vupper3 = Int[]
    for i in 1:lastindex(m_temp)-1
        if m_temp[i] > m_temp[i+1]
            m_temp = deleteat!(m_temp, i+1)
            v_temp = deleteat!(v_temp, i+1)
            c_temp = deleteat!(c_temp, i+1)
            return(m_temp, c_temp, v_temp)
            break
        end
    end
    return(m_temp, c_temp, v_temp)
end

function second_upper_envelope(m_temp, c_temp, v_temp) # Repeat SUE until vectors are clear
    length_original = length(m_temp)
    second_upper_envelope_step(m_temp,c_temp ,v_temp)
    length_new = length(m_temp)
    while length_original > length_new
        length_original = length(m_temp)
        second_upper_envelope_step(m_temp, c_temp, v_temp)     
        length_new = length(m_temp)
    end
    return(m_temp, c_temp, v_temp)
end

function second_upper_envelope_step_pur(m_temp, c_temp, h_temp, v_temp) # Second upper envelope
    for i in 1:lastindex(m_temp)-1
        if m_temp[i] > m_temp[i+1]
            m_temp = deleteat!(m_temp, i+1)
            v_temp = deleteat!(v_temp, i+1)
            c_temp = deleteat!(c_temp, i+1)
            h_temp = deleteat!(h_temp, i+1)
            return(m_temp, c_temp, h_temp, v_temp)
            break
        end
    end
    return(m_temp, c_temp, h_temp, v_temp)
end

function second_upper_envelope_pur(m_temp, c_temp, h_temp, v_temp) # Repeat SUE until vectors are clear
    length_original = length(m_temp)
    second_upper_envelope_step_pur(m_temp,c_temp, h_temp, v_temp)
    length_new = length(m_temp)
    while length_original > length_new
        length_original = length(m_temp)
        second_upper_envelope_step_pur(m_temp, c_temp, h_temp, v_temp)     
        length_new = length(m_temp)
    end
    return(m_temp, c_temp, h_temp, v_temp)
end

function dcegm_s(grid, y, o, h, τ, k, θ) # Sellers

    ψ_o = ψ_r

    eopa_range = grid 

    c_temp = zeros(length(grid))
    m_temp = zeros(length(grid))
    v_temp = zeros(length(grid))
    rh = zeros(length(grid))
    EV = zeros(length(grid))

    c0 = zeros(length(grid))
    h0 = zeros(length(grid))
    v0 = zeros(length(grid))

    for (i,r) in enumerate(eopa_range) 
        
        # Initialize rh and EV value for ith eopa 
        rh_i = 0 
        EV_i = 0

        # Iterate over each state of nature to find future levels of assets
        if T-τ < N
            for ys = 1:lastindex(y_shocks_perm[:,1])
                for rh = 1:lastindex(rh_backw[1,:])

                    pt1 = pt[T-τ] * (1 + rh_backw[rh+2, 1])

                    Yt1 = Y[y, T-τ]*(1.0 + y_shocks_perm[ys,1])

                    Lt1 = r * (1+r_free) + h_bar[h] * ω[o] * pt1      

                    if Lt1 < 0 
                        # Aggregate RHS
                        rh_i += 0
                        # Aggregate EV
                        EV_i += -10e1/ (y_shocks_perm[ys,2]*rh_backw[rh+2,2]) 
                    else
                        cplus_r = C_E_r[T-τ+1, : , :] 
                        citp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), cplus_r, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                        Ct1_r = citp_r(Yt1, Lt1)

                        vplus_r = V_E_r[T-τ+1, :, :]
                        vitp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), vplus_r, Gridded(Linear())), Interpolations.Flat())
                        Vt1_r = vitp_r(Yt1, Lt1)

                        # Now you can calculate the RHS for this specific state of nature 
                        rhs_r = (1 -π1[T-τ]) * β * (1 + r_free)  * q(pt1) * (Ct1_r * q(pt1)) ^(-γ) + 
                        π1[T-τ]*θ* β * (1 + r_free) * (Lt1 + k) ^(-γ)
                        
                        # Aggregate RHS
                        rh_i += rhs_r/(y_shocks_perm[ys,2]*rh_backw[rh+2,2])

                        # Aggregate EV
                        EV_i += (((1 -π1[T-τ]) * β *  Vt1_r) +π1[T-τ]*θ* β * v_w(Lt1, 0, 0, pt1, k))/ (y_shocks_perm[ys,2]*rh_backw[rh+2,2])
                    end
                end
            end
        elseif T-τ >= N && T-τ < T
            for rh = 1:lastindex(rh_backw[1,:])

                pt1 = pt[T-τ] * (1 + rh_backw[rh+2, 1])

                if T-τ == N
                    Yt1 =  Y[y, T-τ] * ϕ 
                else 
                    Yt1 =  Y[y, T-τ]
                end 

                Lt1 = r * (1+r_free) +  h_bar[h] * ω[o] * pt1

                if Lt1 < 0 
                    # Aggregate RHS 
                    rh_i += 0
                    # Aggregate EV
                    EV_i += -10e1/(rh_backw[rh+2,2])
                else
                cplus_r = C_E_r[T-τ+1, : , :] 
                citp_r = extrapolate(interpolate((Y[:, T-τ+1] , grid), cplus_r, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                Ct1_r = citp_r(Yt1, Lt1)

                vplus_r = V_E_r[T-τ+1, :, :]
                vitp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), vplus_r, Gridded(Linear())), Interpolations.Flat())
                Vt1_r = vitp_r(Yt1, Lt1)

                # Now you can calculate the RHS for this specific state of nature 
                rhs_r = (1 -π1[T-τ]) * β * (1 + r_free)  * q(pt1) * (Ct1_r * q(pt1)) ^(-γ) + 
                   π1[T-τ]*θ* β * (1 + r_free)   * (Lt1 + k) ^(-γ)
                
                # Aggregate RHS
                rh_i += rhs_r/(rh_backw[rh+2,2])
            
                # Aggregate EV
                EV_i += ((1 -π1[T-τ]) * β * Vt1_r +π1[T-τ] * θ * β * v_w(Lt1, 0, 0, pt1, k))/(rh_backw[rh+2,2])
                end
            end
        end 

        rh[i] = rh_i
        EV[i] = EV_i

        expon = 1/((1-η) * (1-γ) - 1)
        if rh[i] > 0 
            c_temp[i]  = (h_bar[h] ^ (-η * (1-γ)) * (rh[i] * inv_eta)) ^ (expon) # Endogenous level of consumption
            m_temp[i] = r + c_temp[i] + h_bar[h] *  pt[T-τ] * ((1-ω[o])*ψ_o) + h_bar[h] * pt[T-τ] * sell_costs  # Endogenous level of resources
            v_temp[i] =  u_c(c_temp[i], h_bar[h], pt[T-τ]) + EV[i] 
        else 
            c_temp[i] = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] = -10e1
        end
        if m_temp[i] < 0 || c_temp[i] < 0
            c_temp[i] = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] = -10e1
        end
    end

    unique!(c_temp)
    unique!(m_temp)
    unique!(v_temp)

    second_upper_envelope(m_temp, c_temp, v_temp) 

    # Extrapolate consumption based on common grid for liquid assets
    c_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [0, c_temp...], Gridded(Linear())), Interpolations.Flat()) # Recover consumption levels from endogenous grid
    v_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [-10e1, v_temp...], Gridded(Linear())), Interpolations.Flat()) 
    
    for (i,r) in enumerate(Y[y, T-τ] .+  grid) 
        c0[i] = c_extrapol(r)
        v0[i] = v_extrapol(r)
        if v0[i] == NaN || v0[i] == -Inf
            v0[i] = -10e1
        end
    end
    return(c0, v0)
end

function dcegm_o(grid, y, o, h, τ, k, θ) # Owners

    ψ_o = ψ_r

    eopa_range = grid 

    c_temp = zeros(length(grid))
    m_temp = zeros(length(grid))
    v_temp = zeros(length(grid))
    rh = zeros(length(grid))
    EV = zeros(length(grid))

    c0 = zeros(length(grid))
    v0 = zeros(length(grid))

    for (i,r) in enumerate(eopa_range) 
        #Initialize rh and EV value for ith eopa
        rh_i = 0 
        EV_i  = 0 
        # Iterate over each state of nature to find future levels of assets
        if T-τ < N
            for ys = 1:lastindex(y_shocks_perm[:,1])
                #for rt = 1:lastindex(y_shocks_temp[:,1])
                    for rh = 1:lastindex(rh_backw[1,:])
                        Yt1 = Y[y, T-τ]*(1.0 + y_shocks_perm[ys,1])
                        Lt1 = r * (1+r_free)
                        #Lt1 = Lt1*(1.0 + y_shocks_temp[rt,1])
                        pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                        ω1 = (pt1/pt[T-τ] - (1-ω[o]))/(pt1/pt[T-τ])
                        if Yt1 + Lt1 - h_bar[h] * pt1 * ((1-ω1)*ψ_o)  > 0  && Lt1 + h_bar[h] * ω1 * pt1 >= 0

                             # If we remain owners
                            cplus_p = C_E_o[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                            citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_p, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                            Ct1_p = citp_p(Yt1, Lt1, ω1)

                            vplus_p = V_E_o[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                            vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_p, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                            Vt1_p = vitp_p(Yt1, Lt1, ω1)

                            rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                            π1[T-τ]*θ* β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)

                            # If we sell
                            Lt1s = Lt1 -  h_bar[h]*pt1*sell_costs
                            cplus_s = C_E_s[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                            citp_s = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_s, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                            Ct1_s = citp_s(Yt1, Lt1s, ω1)

                            vplus_s = V_E_s[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                            vitp_s = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_s, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                            Vt1_s = vitp_s(Yt1, Lt1s, ω1)
                            rhs_s = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_s)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                            π1[T-τ]*θ* β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)

                            # Aggregate RHS
                            rh_i += maximum(last, (rhs_p, rhs_s))/ (rh_backw[rh,2]*y_shocks_perm[ys,2])
                            
                             # Aggregate EV
                            EV_i += ((1 -π1[T-τ]) * β *   maximum(last, (Vt1_p, Vt1_s)) +π1[T-τ] * θ * β * v_w(Lt1, ω1, h, τ, k))/(rh_backw[rh,2]*y_shocks_perm[ys,2])
                            
                        else
                            rh_i += 0
                            EV_i += -10e1/ (rh_backw[rh,2]*y_shocks_perm[ys,2])
                        end
                    end 
                #end
            end
        elseif  T-τ >= N && T-τ < T-1
            for rh = 1:lastindex(rh_backw[1,:])
                if T-τ == N
                    Yt1 =  Y[y, T-τ] * ϕ 
                else 
                    Yt1 =  Y[y, T-τ]
                end
                Lt1 = r * (1+r_free)
                pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                ω1 = (pt1/pt[T-τ] - (1-ω[o]))/(pt1/pt[T-τ])
                if Yt1 + Lt1 - h_bar[h] * pt1 * ((1-ω1)*ψ_o)  > 0  && Lt1 + h_bar[h] * ω1 * pt1 >= 0

                    # If we remain owners
                    cplus_p = C_E_o[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                    citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_p, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                    Ct1_p = citp_p(Yt1, Lt1, ω1)

                    vplus_p = V_E_o[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                    vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_p, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                    Vt1_p = vitp_p(Yt1, Lt1, ω1)

                    rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                    π1[T-τ]*θ* β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)

                    # If we sell
                    Lt1s = Lt1 -  h_bar[h]*pt1*sell_costs
                    cplus_s = C_E_s[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                    citp_s = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_s, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                    Ct1_s = citp_s(Yt1, Lt1s, ω1)

                    vplus_s = V_E_s[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                    vitp_s = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_s, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                    
                    vplus_s = V_E_s[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                    vitp_s = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_s, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                    Vt1_s = vitp_s(Yt1, Lt1s, ω1)
                    rhs_s = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_s)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                    π1[T-τ]*θ* β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)

                    # Aggregate RHS
                    rh_i += maximum(last, (rhs_p, rhs_s))/(rh_backw[rh,2])
                    
                    # Aggregate EV
                    EV_i += ((1 -π1[T-τ]) * β *   maximum(last, (Vt1_p, Vt1_s)) +π1[T-τ] * θ * β * v_w(Lt1, ω1, h, τ, k))/(rh_backw[rh,2])
                        
                else
                    rh_i += 0
                    EV_i += -10e1/(rh_backw[rh,2]) 
                end
            end
        elseif T-τ == T-1
            for rh = 1:lastindex(rh_backw[1,:])
                Yt1 = Y[y, T-τ+1]
                Lt1 = r * (1+r_free)
                pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                ω1 = (pt1/pt[T-τ] - (1-ω[o]))/(pt1/pt[T-τ])
                if Yt1 + Lt1 - h_bar[h] * pt1 * ((1-ω1)*ψ_o)  > 0 && Lt1 + h_bar[h] * ω1 * pt1 >= 0

                    cplus_p = C_E_o[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                    citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_p, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                    Ct1_p = citp_p(Yt1, Lt1, ω1)

                    vplus_p = V_E_o[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                    vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_p, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                    Vt1_p = vitp_p(Yt1, Lt1, ω1)

                    # Now you can calculate the RHS for this specific state of nature 
                    rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                       π1[T-τ]*θ* β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)
                    
                        # Aggregate RHS
                    rh_i += rhs_p/(rh_backw[rh,2])

                    # Aggregate EV
                    EV_i += ((1 -π1[T-τ]) * β *  Vt1_p +π1[T-τ]*θ* β * v_w(Lt1, ω1, h, pt1, k))/(rh_backw[rh,2])
                else     
                    rh_i += 0
                    EV_i += -10e1/(rh_backw[rh,2]) 
                end
            end
        elseif T-τ == T
            for rh = 1:lastindex(rh_backw[1,:])
                Lt1 = r * (1+r_free)
                pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                ω1 = (pt1/pt[T-τ] - (1-ω[o]))/(pt1/pt[T-τ])

                if Lt1 + h_bar[h] *  ω1 * pt1 >= 0 
                    # Now you can calculate the RHS for this specific state of nature 
                    rhs_p = θ * β * (1 + r_free)  * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)

                    # Aggregate RHS
                    rh_i += rhs_p/(rh_backw[rh,2])

                    # Aggregate EV
                    EV_i += (θ* β * v_w(Lt1, ω1, h, pt1, k))/(rh_backw[rh,2])
                else
                    rh_i += 0
                    EV_i += -10e1/(rh_backw[rh,2]) 
                end
            end
        end

        rh[i] = rh_i
        EV[i] = EV_i
        expon = 1/((1-η) * (1-γ) - 1)
        if rh[i] > 0 
            c_temp[i]  = (h_bar[h] ^ (-η * (1-γ)) * (rh[i] * inv_eta)) ^ (expon) # Endogenous level of consumption
            m_temp[i] = r + c_temp[i] + h_bar[h] * pt[T-τ] * ((1-ω[o])*ψ_o) # Endogenous level of resources
            v_temp[i] =  u_c(c_temp[i], h_bar[h], pt[T-τ]) + EV[i] 
        else 
            c_temp[i] = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] =  -10e1
        end
        if m_temp[i] < 0 || c_temp[i] < 0  || v_temp[i] < -10e1
            c_temp[i] = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] = -10e1
        end
    end

    m_temp = [x for x in m_temp if x > 0.0001]
    c_temp = [x for x in c_temp if x > 0.0001]
    v_temp = [x for x in v_temp if x > -10e1]
        
    unique!(c_temp)
    unique!(m_temp)
    unique!(v_temp)

    second_upper_envelope(m_temp, c_temp, v_temp)

    # Extrapolate consumption based on common grid for liquid assets
    c_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [0, c_temp...], Gridded(Linear())), Interpolations.Flat()) # Recover consumption levels from endogenous grid
    v_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [-10e1, v_temp...], Gridded(Linear())), Interpolations.Flat()) 

    for (i,r) in enumerate(Y[y, T-τ] .+  grid) 
        c0[i] = c_extrapol(r)
        v0[i] = v_extrapol(r)
        if v0[i] == NaN || v0[i] == -Inf
            v0[i] = -10e1
        end
    end
    
    return(c0, v0)    
end

function dcegm_r(grid, y, τ, k, θ) # Renters

    eopa_range = grid 

    c_temp = zeros(length(grid))
    v_temp = zeros(length(grid))
    m_temp = zeros(length(grid))
    rh = zeros(length(grid))
    EV = zeros(length(grid))
    c0 = zeros(length(grid))
    h0 = zeros(length(grid))
    v0 = zeros(length(grid))

    for (i,r) in enumerate(eopa_range) # For each EOPA we  need the RHS for both cases (renter, purchaser)
        # Initialize rh and EV value for ith eopa
        rh_i = 0 
        EV_i  = 0 

        # Iterate over each state of nature to find future levels of assets
        if T-τ < N
            for ys = 1:lastindex(y_shocks_perm[:, 1])
                #for rt = 1:lastindex(y_shocks_temp[:,1])
                    Yt1 = Y[y, T-τ]*(1.0 + y_shocks_perm[ys,1])
                    Lt1 = r * (1+r_free)

                    if  Lt1 < 0
                        EV_i += -10e1/ (y_shocks_perm[ys,2]) 
                        rh_i += 0

                    else
                        for rh = 1:lastindex(rh_backw[1,:])
                            pt1 = pt[T-τ] * (1 + rh_backw[rh+2, 1])

                            # If we stay renters:
                            cplus_r = C_E_r[T-τ+1, :, :]
                            citp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), cplus_r, Gridded(Linear())), Interpolations.Flat())
                            Ct1_r = citp_r(Yt1, Lt1)

                            vplus_r = V_E_r[T-τ+1, :, :]
                            vitp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), vplus_r, Gridded(Linear())), Interpolations.Flat())
                            Vt1_r = vitp_r(Yt1, Lt1)
                        
                            rhs_r = (1 -π1[T-τ]) * β * (1 + r_free)  * q(pt1) * (Ct1_r * q(pt1)) ^(-γ) + 
                                   π1[T-τ]*θ* β * (1 + r_free)   * (Lt1 + k) ^(-γ)
                            
                            # If we purchase
                            cplus_p = C_E_p[T-τ+1, :, :] # Get grid of future consumption for all st and wrh levels 
                            vplus_p = V_E_p[T-τ+1, :, :]
                            hplus_p = H_E_p[T-τ+1, :, :]
                            citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid), cplus_p, Gridded(Linear())), Interpolations.Flat())
                            vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid), vplus_p, Gridded(Linear())), Interpolations.Flat())
                            hitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid), hplus_p, Gridded(Constant())), Interpolations.Flat())
                            Ht1_p = hitp_p(Yt1, Lt1)
                            Ct1_p = citp_p(Yt1, Lt1)
                            if Yt1 + Lt1 - Ct1_p - Ht1_p*pt1*((ψ_r*(1-δ)+δ)) - Ht1_p * pt1 * sell_costs1 < 0
                                Vt1_p = -10e1
                            else
                                Vt1_p = vitp_p(Yt1, Lt1)
                            end
                            if Vt1_p <= -10e1 || Ht1_p < 30
                                rhs_p = 0  
                            else 
                                rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * Ht1_p^(η*(1-γ))  +
                               π1[T-τ]*θ* β * (1 + r_free)   * (Lt1 + k) ^(-γ)
                            end
                            
                            EV_i += ((1 -π1[T-τ]) * β * maximum(last, (Vt1_p, Vt1_r)) +π1[T-τ] * θ * β * v_w(Lt1, 0, 0, pt1, k))/(rh_backw[rh+2,2]*y_shocks_perm[ys,2])

                            # Aggregate RHS
                            rh_i += maximum(last, (rhs_p, rhs_r))/(rh_backw[rh+2,2]*y_shocks_perm[ys,2])
                        end
                    end
                #end
            end
        elseif  T-τ >= N && T-τ < T-1
            if T-τ == N
                    Yt1 =  Y[y, T-τ] * ϕ 
            else 
                    Yt1 =  Y[y, T-τ]
            end 

            Lt1 = r * (1+r_free) # Future resources 
            
            if  Lt1 < 0
                EV_i += -10e1
                rh_i += 0

            else
                for rh = 1:lastindex(rh_backw[1,:])
                    pt1 = pt[T-τ] * (1 + rh_backw[rh+2, 1])

                    # If we stay renters:
                    cplus_r = C_E_r[T-τ+1, :, :]
                    citp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), cplus_r, Gridded(Linear())), Interpolations.Flat())
                    Ct1_r = citp_r(Yt1, Lt1)

                    vplus_r = V_E_r[T-τ+1, :, :]
                    vitp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), vplus_r, Gridded(Linear())), Interpolations.Flat())
                    Vt1_r = vitp_r(Yt1, Lt1)
                
                    rhs_r = (1 -π1[T-τ]) * β * (1 + r_free)  * q(pt1) * (Ct1_r * q(pt1)) ^(-γ) + 
                           π1[T-τ]*θ* β * (1 + r_free)   * (Lt1 + k) ^(-γ)
                    
                    # If we purchase
                    cplus_p = C_E_p[T-τ+1, :, :] # Get grid of future consumption for all st and wrh levels 
                    vplus_p = V_E_p[T-τ+1, :, :]
                    hplus_p = H_E_p[T-τ+1, :, :]
                    citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid), cplus_p, Gridded(Linear())), Interpolations.Flat())
                    vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid), vplus_p, Gridded(Linear())), Interpolations.Flat())
                    hitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid), hplus_p, Gridded(Constant())), Interpolations.Flat())
                    Ht1_p = hitp_p(Yt1, Lt1)
                    Ct1_p = citp_p(Yt1, Lt1)
                    if Yt1 + Lt1 - Ct1_p - Ht1_p*pt1*((ψ_r*(1-δ)+δ)) - Ht1_p * pt1 * sell_costs1 < 0
                        Vt1_p = -10e1
                    else
                        Vt1_p = vitp_p(Yt1, Lt1)
                    end
                    if Vt1_p <= -10e1 || Ht1_p < 30
                        rhs_p = 0  
                    else 
                        rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * Ht1_p^(η*(1-γ))  +
                       π1[T-τ]*θ* β * (1 + r_free) * (Lt1 + k) ^(-γ)
                    end
                    
                    EV_i += ((1 -π1[T-τ]) * β * maximum(last, (Vt1_p, Vt1_r)) +π1[T-τ] * θ * β * v_w(Lt1, 0, 0, pt1, k))/(rh_backw[rh+2,2])

                    # Aggregate RHS
                    rh_i += maximum(last, (rhs_p, rhs_r))/(rh_backw[rh+2,2])
                end
            end
        elseif T-τ == T -1 
            Yt1 =  Y[y, T-τ]
            Lt1 = r * (1+r_free)
            if  Lt1 < 0
                # Aggregate RHS
                rh_i += 0
                # Aggregate EV
                EV_i += -10e1
            else
                for rh = 1:lastindex(rh_backw[1,:])
                    pt1 = pt[T-τ] * (1 + rh_backw[rh+2, 1])
                    cplus_r = C_E_r[T-τ+1, : , :] 
                    citp_r = extrapolate(interpolate((Y[:, T-τ+1] , grid), cplus_r, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                    Ct1_r = citp_r(Yt1, Lt1)

                    vplus_r = V_E_r[T-τ+1, :, :]
                    vitp_r = extrapolate(interpolate((Y[:, T-τ+1], grid), vplus_r, Gridded(Linear())), Interpolations.Flat())
                    Vt1_r = vitp_r(Yt1, Lt1)

                    # Now you can calculate the RHS for this specific state of nature 
                    rhs_r = (1 -π1[T-τ]) * β * (1 + r_free)  * q(pt1) * (Ct1_r * q(pt1)) ^(-γ) + 
                       π1[T-τ]*θ* β * (1 + r_free)   * (Lt1 + k) ^(-γ)
                    
                    # Aggregate RHS
                    rh_i += rhs_r/(rh_backw[rh+2,2])
                
                    # Aggregate EV
                    EV_i += ((1 -π1[T-τ]) * β * Vt1_r +π1[T-τ] * θ * β * v_w(Lt1, 0, 0, pt1, k))/(rh_backw[rh+2,2])
                end
            end
        elseif T-τ ==  T
            Lt1 = r * (1+r_free)
            if Lt1 < 0
                # Aggregate RHS
                rh_i += 0
                # Aggregate EV
                EV_i += -10e1
            else                        
                # Now you can calculate the RHS for this specific state of nature 
                rhs_r =π1[T-τ]*θ* β * (1 + r_free) * (Lt1 + k) ^(-γ)
                
                # Aggregate RHS
                rh_i += rhs_r

                # Aggregate EV
                EV_i += (π1[T-τ]*θ* β * v_w(Lt1, 0, 0, 0, k))
            end
        end 

        rh[i] = rh_i
        EV[i] = EV_i
        
        if  rh[i] > 0 
            c_temp[i] = ((1/q(pt[T-τ])^(1-γ)) * rh[i]) ^ (-1/γ) # Apply inverse euler 
            m_temp[i] = r + c_temp[i] + h_c(c_temp[i], 0, τ) * ψ_r * pt[T-τ] # Endogenous level of resources
            v_temp[i] = u_c(c_temp[i], 0, pt[T-τ]) + EV[i]
        else
            c_temp[i]  = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] = -10e1 
        end
        if m_temp[i] < 0 || c_temp[i] < 0  || v_temp[i] < -10e1
            c_temp[i] = 1/10000
            m_temp[i] = 1/10000
            v_temp[i] = -10e1
        end
    end

    unique!(c_temp)
    unique!(m_temp)
    unique!(v_temp)

    second_upper_envelope(m_temp, c_temp, v_temp)

    # Extrapolate consumption based on common grid for liquid assets
    c_extrapol = extrapolate(
       interpolate(([0, m_temp...],), [0, c_temp...], Gridded(Linear())), Interpolations.Flat()) # Recover consumption levels from endogenous grid
    
    v_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [-10e1, v_temp...], Gridded(Linear())), Interpolations.Flat()) 

    for (i,r) in enumerate(Y[y, T-τ] .+ grid)
        c0[i] = c_extrapol(r)
        h0[i] = h_c(c0[i], 0, τ)
        v0[i] = v_extrapol(r) 
    end
    return(c0, h0, v0)
end

function dcegm_p(grid, y, τ, k, θ) # Purchasers

    ψ_o = ψ_r

    eopa_range = grid
    
    c_temp = zeros(length(grid))
    h_temp = zeros(length(grid))
    m_temp = zeros(length(grid))
    v_temp = zeros(length(grid))

    rht = zeros(length(grid))
    
    c0 = zeros(length(grid))
    v0 = zeros(length(grid))
    h0 = zeros(length(grid))

    for (i,r) in enumerate(eopa_range)
        c_vec = zeros(length(h_bar))
        m_vec = zeros(length(h_bar))
        v_vec = zeros(length(h_bar))
        rhth = zeros(length(h_bar))
        for h = 1:lastindex(h_bar) # Iterate over discrete values of real estate
            r_net = r
            if r_net > 0 
                #initialize rh and EV value for ith eopa
                rh_i = 0 
                EV_i  = 0 
                
                # Iterate over each state of nature to find future levels of assets
                if T-τ < N
                    for ys = 1:lastindex(y_shocks_perm[:,1]) 
                        #for rt = 1:lastindex(y_shocks_temp[:,1])
                            for rh = 1:lastindex(rh_backw[1,:])
                                Yt1 = Y[y, T-τ]*(1.0 + y_shocks_perm[ys,1])
                                Lt1 = r_net * (1+r_free)
                                pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                                ω1 = (pt1/pt[T-τ] - (1-δ))/(pt1/pt[T-τ])
                                if  Yt1 + Lt1 - h_bar[h] * pt1 * ((1-ω1)*ψ_o) > 0 && Lt1 + h_bar[h] * ω1 * pt1 >= 0 
                                    cplus_p = C_E_o[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                                    citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_p, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                                    Ct1_p = citp_p(Yt1, Lt1, ω1)

                                    vplus_p = V_E_o[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                                    vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_p, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                                    Vt1_p = vitp_p(Yt1, Lt1, ω1)

                                    # Now you can calculate the RHS for this specific state of nature 
                                    rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                                   π1[T-τ]*θ* β * (1 + r_free) * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)
                                    
                                    # Aggregate RHS
                                    rh_i += rhs_p/(rh_backw[rh,2]*y_shocks_perm[ys,2])

                                    # Aggregate EV
                                    EV_i += ((1 -π1[T-τ]) * β *  Vt1_p +π1[T-τ] * θ * β * v_w(Lt1, ω1, h, pt1, k))/(rh_backw[rh,2]*y_shocks_perm[ys,2]) 
                                else 
                                    rh_i += 0
                                    EV_i += -10e1/ (rh_backw[rh,2]*y_shocks_perm[ys,2]) 
                                end
                            end
                        #end
                    end
                elseif T-τ >= N
                    for rh = 1:lastindex(rh_backw[1,:])

                        if T-τ == N
                            Yt1 =  Y[y, T-τ] * ϕ 
                        else 
                            Yt1 =  Y[y, T-τ]
                        end 

                        Lt1 = r_net * (1+r_free)
                        pt1 = pt[T-τ]*(1 + rh_backw[rh, 1]) 
                        ω1 = (pt1/pt[T-τ] - (1-δ))/(pt1/pt[T-τ])
                        if  Yt1 + Lt1 - h_bar[h] * pt1 * ((1-ω1)*ψ_o) > 0 && Lt1 + h_bar[h] * ω1 * pt1 >= 0
                            cplus_p = C_E_o[T-τ+1, :, :, :, h] # Get grid of future consumption for all st and wrh levels 
                            citp_p = extrapolate(interpolate((Y[:, T-τ+1], grid, ω[:],), cplus_p, Gridded(Linear())), Interpolations.Flat()) # Consumption function based on future utility levels (as you use the inverse marginal utility) 
                            Ct1_p = citp_p(Yt1, Lt1, ω1)

                            vplus_p = V_E_o[T-τ+1, :, :, :, h]  # Get grid of future utility for all st and wrh levels 
                            vitp_p = extrapolate(interpolate((Y[:, T-τ+1] , grid, ω[:],), vplus_p, Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
                            Vt1_p = vitp_p(Yt1, Lt1, ω1)

                            # Now you can calculate the RHS for this specific state of nature 
                            rhs_p = (1 -π1[T-τ]) * β * (1 + r_free) * (1-η) * (Ct1_p)^((1-η)*(1-γ)-1) * h_bar[h]^(η*(1-γ))  +
                           π1[T-τ]*θ* β * (1 + r_free) * (Lt1 + h_bar[h] * ω1 * pt1 + k) ^(-γ)
                            
                            # Aggregate RHS
                            rh_i += rhs_p/(rh_backw[rh,2])

                            # Aggregate EV
                            EV_i += ((1 -π1[T-τ]) * β *  Vt1_p +π1[T-τ] * θ * β * v_w(Lt1, ω1, h, pt1, k))/(rh_backw[rh,2])
                        else 
                            rh_i += 0
                            EV_i += -10e1/(rh_backw[rh,2]) 
                        end
                    end
                end
                if rh_i <= 0 
                    rh_i = 0 
                    EV_i = -10e1
                    c_vec[h] = 0 
                    m_vec[h] = 0
                    v_vec[h] = -10e1
                    rhth[h] = 0 
                else
                    expon = 1/((1-η) * (1-γ) - 1)
                    c_vec[h]  = (h_bar[h] ^ (-η * (1-γ)) * (rh_i * inv_eta)) ^ (expon)  
                    m_vec[h] = r_net + c_vec[h] + h_bar[h] * ((1-δ)*ψ_o+δ) * pt[T-τ] + h_bar[h] * pt[T-τ] * sell_costs1 # Endogenous level of resources
                    v_vec[h] = u_c(c_vec[h], h_bar[h], pt[T-τ]) + EV_i
                    rhth[h] = rh_i
                end
            else
                v_vec[h] = -10e1
            end
        end
        
        mi = findmax(v_vec)[2]
        c_temp[i]  = c_vec[mi]
        rht[i] = rhth[mi]
        m_temp[i] =  m_vec[mi]
        h_temp[i] =  h_bar[mi]
        v_temp[i] =  v_vec[mi]

    end
    
    matr = [m_temp c_temp h_temp v_temp]

    # Extract the first column
    first_column = matr[:, 1]

    # Find unique values in the first column and their indices
    unique_values =  unique(i -> first_column[i], 1:length(first_column))
    # Use the indices to get the rows with unique first column values
    matr = matr[unique_values, :]

    m_temp = matr[:, 1]
    c_temp = matr[:, 2]
    h_temp = matr[:, 3]
    v_temp = matr[:, 4]

   second_upper_envelope_pur(m_temp, c_temp, h_temp, v_temp)

    # Extrapolate consumption based on common grid for liquid assets
    c_extrapol = extrapolate(
       interpolate(([0, m_temp...],), [0, c_temp...], Gridded(Linear())), Interpolations.Flat()) # Recover consumption levels from endogenous grid
    h_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [0, h_temp...], Gridded(Constant())), Interpolations.Flat()) 
    v_extrapol = extrapolate(
        interpolate(([0, m_temp...],), [-10e1, v_temp...], Gridded(Linear())), Interpolations.Flat()) 
    
    for (i,r) in enumerate(Y[y, T-τ] .+ grid)
        c0[i] = c_extrapol(r)
        h0[i] = h_extrapol(r)
        v0[i] = v_extrapol(r) 

       if c0[i] > r
            c0[i] = r
            h0[i] = 0
        end
   end
    return(c0, h0, v0)
end

###########################################################################################################
function prob_d_f(VR, VP) # Purchase decision
    if VR>=VP 
        p0 = 1
    else
        p0 = 0       
    end
    return(p0)
end

function closest_value_s(t, val) # Refinement functions for defaulting
    dxfirst = val - S[1]
    dxindex = 1
    for (i,r) in enumerate(S[:])
        dx = val - r 
        if dx > 0 && dx < dxfirst
            dxindex = i
        end
    end
    return(S[dxindex])
end 

function closest_value_y(t, val)
    dxfirst = val - Y[1,t]
    dxindex = 1
    for (i,r) in enumerate(Y[:,t])
        dx = val - r 
        if dx > 0 && dx < dxfirst
            dxindex = i
        end
    end
    return(Y[dxindex,t])
end 

begin # Load data for housing prices
    vct = 1:T

    ptmat = vcat(vct', pt')'
    
    ptmat = DataFrame(ptmat, :auto)
    rename!(ptmat, "x1" => "time", "x2" => "pt")
    ptmat.time = Int.(ptmat.time)
end

function df_creation(agent_n, seed_n, shock_perm, shock_temp, save) # Main df
    A = agent_n
    ψ_o = ψ_r
    ####################################################################################################################
    ### Create matrices
    YT = zeros(Float64,(A, T)) # create matrix of income profiles
    CT = zeros(Float64,(A, T)) # create matrix of housing profiles
    HT = zeros(Float64,(A, T)) # create matrix of housing profiles
    ST = zeros(Float64,(A, T)) # create matrix of liquid assets profiles
    LT = zeros(Float64,(A, T)) # create matrix of cash at hand profiles
    VT = zeros(Float64,(A, T)) # create matrix of utility profiles
    PR = zeros(Float64,(A, T)) # create matrix of purchasing probability
    OT = zeros(Float64,(A, T)) # create matrix of mortgages
    DT = zeros(Float64,(A, T)) # create matrix for ownership
    DP = zeros(Float64,(A, T)) # create matrix for purchase
    DS = zeros(Float64,(A, T)) # create matrix for selling
    PP = zeros(Float64,(A, T))  # create matrix for price at time of purchase
    SH = zeros(Float64,(A, T))  # create matrix for time of the shock
    DF = zeros(Float64,(A, T))  # create matrix for default
    ####################################################################################################################
    # Life-Cycle simulation
    nr = 0 

    for i in 1:A 
        SH[i, :] .= Int(vect_sh[i])
    end

    for t in 1:T
        ## Interpolations
        #### NON PARTICIPATING
        # Renters
        cr_extp_t = extrapolate(interpolate((Y[:, t],S[:],), C_E_r[t, :, :], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
        hr_extp_t = extrapolate(interpolate((Y[:, t],S[:],), H_E_r[t, :, :], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
        vr_extp_t = extrapolate(interpolate((Y[:, t],S[:],), V_E_r[t, :, :], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation

        # purchasers
        cp_extp_t = extrapolate(interpolate((Y[:, t],S[:],), C_E_p[t, :, :], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation
        hp_extp_t = extrapolate(interpolate((Y[:, t],S[:],), H_E_p[t, :, :], Gridded(Constant())), Interpolations.Flat())  # Two-dimensional interpolation
        vp_extp_t = extrapolate(interpolate((Y[:, t],S[:],), V_E_p[t, :, :], Gridded(Linear())), Interpolations.Flat())  # Two-dimensional interpolation

        # Owners
        co_extp_t = extrapolate(interpolate((Y[:, t],S[:], ω[:], h_bar,), C_E_o[t, :, :, :, :], Gridded(Linear())), Interpolations.Flat()) # Four-dimensional interpolation
        vo_extp_t = extrapolate(interpolate((Y[:, t],S[:], ω[:], h_bar,), V_E_o[t, :, :, :, :], Gridded(Linear())), Interpolations.Flat())  # Four-dimensional interpolation

        # Sellers
        cs_extp_t = extrapolate(interpolate((Y[:, t],S[:], ω[:], h_bar,), C_E_s[t, :, :, :, :], Gridded(Linear())), Interpolations.Flat()) # Four-dimensional interpolation
        vs_extp_t = extrapolate(interpolate((Y[:, t],S[:], ω[:], h_bar,), V_E_s[t, :, :, :, :], Gridded(Linear())), Interpolations.Flat())  # Four-dimensional interpolation
        
        for i in 1:A
            if t == 1
                YT[i, 1] = YT_init[i]
            else
                if t <= N
                    if shock_perm == 0 
                        YT[i, t] = YT[i, t-1]*rn[i*t]*ypsilon[t]
                    elseif  shock_perm == 1 && SH[i, t] != float(t)
                        YT[i, t] = YT[i, t-1]*rn[i*t]*ypsilon[t]
                    elseif  shock_perm == 1 && SH[i, t] == float(t)
                        YT[i, t] = YT[i, t-1]*0.7
                    end
                else
                    YT[i, t] = YT[i, N]* ϕ
                end

                LT[i, t] =  ST[i, t-1]*(1+r_free)
                if LT[i, t] < 1/10000
                    LT[i, t] = 0 
                end

                if DT[i, t] == 1 && DP[i, t] != 1
                    OT[i, t] = (pt[t]/pt[t-1] - (1-OT[i, t-1]))/(pt[t]/pt[t-1])
                end

                if DS[i, t-1] == 1 && DF[i, t-1] == 0
                    LT[i, t] =  LT[i, t] + OT[i, t-1]*HT[i, t-1]*pt[t]
                end 
            end

            if DT[i, t] == 0 
                if  vp_extp_t(YT[i, t], LT[i, t]) > vr_extp_t(YT[i, t], LT[i, t]) &&
                    hp_extp_t(YT[i, t], LT[i, t]) >=30  &&
                    YT[i, t] + LT[i, t] - cp_extp_t(YT[i, t], LT[i, t]) - hp_extp_t(YT[i, t], LT[i, t])*pt[t]*((ψ_r*(1-δ)+δ+sell_costs1)) > 0

                    nr += 1

                    h_pur = hp_extp_t(YT[i, t], LT[i, t])
                    c_pur = cp_extp_t(YT[i, t], LT[i, t])

                    CT[i, t] =  c_pur
                    HT[i, t] =  h_pur

                    VT[i, t] =  vp_extp_t(YT[i, t], LT[i, t])       
                    ST[i, t] =  YT[i, t] + LT[i, t] -  CT[i, t] - HT[i, t]*pt[t]*((ψ_r*(1-δ)+δ+sell_costs1))

                    if ST[i, t] < 1/10000  && ST[i, t] > -1/10000  # Round small values
                        ST[i, t] = 0 
                    end 

                    OT[i, t] = δ 
                    DT[i, t+1] = 1 
                    DP[i, t]  = 1
                    DT[i, t] = 1
                    PP[i, t:T] .= pt[t]

                else
                    CT[i, t] =  cr_extp_t(YT[i, t], LT[i, t])
                    HT[i, t] =  hr_extp_t(YT[i, t], LT[i, t])
                    VT[i, t] =  vr_extp_t(YT[i, t], LT[i, t])       
                    ST[i, t] =  YT[i, t] + LT[i, t] -  CT[i, t] - HT[i, t]*pt[t]*ψ_r 

                    if ST[i, t] < 1/10000  && ST[i, t] > -1/10000 
                        ST[i, t] = 0 
                    end 
                end

                if ST[i, t] < 0
                    DF[i, t] = 1
                    CT[i, t] =  cr_extp_t(closest_value_y(t, YT[i, t]), closest_value_s(t, LT[i, t]))
                    HT[i, t] =  hr_extp_t(closest_value_y(t, YT[i, t]), closest_value_s(t, LT[i, t]))
                    VT[i, t]  = vr_extp_t(closest_value_y(t, YT[i, t]), closest_value_s(t, LT[i, t]))
                    ST[i, t] = 0
                end

            else 
                HT[i, t] =  HT[i, t-1]   
                # Do we sell?
                if t < T-1
                    if  vs_extp_t(YT[i, t], LT[i, t], OT[i, t], HT[i, t]) > vo_extp_t(YT[i, t], LT[i, t], OT[i, t], HT[i, t]) 
                        DS[i,t] = 1
                        cf = cs_extp_t
                        vf = vs_extp_t
                        DT[i, t+1] = 0 
                    else
                        cf = co_extp_t
                        vf = vo_extp_t
                        DT[i, t+1] = 1 
                    end 
                else
                    cf = co_extp_t
                    vf = vo_extp_t
                    DT[i, t:T] .= 1 
                end

                CT[i, t] =  cf(YT[i, t], LT[i, t], OT[i, t], HT[i, t])
                VT[i, t] =  vf(YT[i, t], LT[i, t], OT[i, t], HT[i, t])    
                ST[i, t] = YT[i, t] + LT[i, t] -  CT[i, t] - HT[i, t]*pt[t]*((1-OT[i, t])*ψ_o) - HT[i, t]*pt[t]*DS[i,t]*sell_costs
                if ST[i, t] < 1/10000 && ST[i, t] > -1/10000 
                    ST[i, t] = 0 
                end 

                if ST[i, t] < 0
                    DF[i, t] = 1
                    CT[i, t] =  YT[i, t] + LT[i, t] - HT[i, t]*pt[t]*((1-OT[i, t])*ψ_o) 
                    HT[i, t] =  HT[i, t-1]                          
                    ST[i, t] = 0 
                    if t < T-1
                        VT[i, t]  =  vs_extp_t(YT[i, t], LT[i, t], OT[i, t], HT[i, t]) 
                        DT[i, t+1] = 0
                        DS[i,t] = 1
                    else 
                        VT[i, t]  =  vf(YT[i, t], LT[i, t], OT[i, t], HT[i, t]) 
                    end
                end
            end
        end
    end
    ####################################################################################################################
    # Stack data and save dataset
    ptmat.time = string.(ptmat.time)

    YT = stack(DataFrame(YT, :auto))
    rename!(YT, "variable" => "time", "value" => "YT")
    YT.time = replace.(YT.time, "x" => "")
    YT.row = rownumber.(eachrow(YT))

    CT = stack(DataFrame(CT, :auto))
    rename!(CT, "variable" => "time", "value" => "CT")
    select!(CT, Not(:time))
    CT.row = rownumber.(eachrow(CT))

    HT = stack(DataFrame(HT, :auto))
    rename!(HT, "variable" => "time", "value" => "HT")
    select!(HT, Not(:time))
    HT.row = rownumber.(eachrow(HT))

    ST = stack(DataFrame(ST, :auto))
    rename!(ST, "variable" => "time", "value" => "ST")
    select!(ST, Not(:time))
    ST.row = rownumber.(eachrow(ST))

    LT = stack(DataFrame(LT, :auto))
    rename!(LT, "variable" => "time", "value" => "LT")
    select!(LT, Not(:time))
    LT.row = rownumber.(eachrow(LT))

    VT = stack(DataFrame(VT, :auto))
    rename!(VT, "variable" => "time", "value" => "VT")
    select!(VT, Not(:time))
    VT.row = rownumber.(eachrow(VT))

    PR = stack(DataFrame(PR, :auto))
    rename!(PR, "variable" => "time", "value" => "PR")
    select!(PR, Not(:time))
    PR.row = rownumber.(eachrow(PR))

    OT = stack(DataFrame(OT, :auto))
    rename!(OT, "variable" => "time", "value" => "OT")
    select!(OT, Not(:time))
    OT.row = rownumber.(eachrow(OT))

    DT = stack(DataFrame(DT, :auto))
    rename!(DT, "variable" => "time", "value" => "DT")
    select!(DT, Not(:time))
    DT.row = rownumber.(eachrow(DT))

    DP = stack(DataFrame(DP, :auto))
    rename!(DP, "variable" => "time", "value" => "DP")
    select!(DP, Not(:time))
    DP.row = rownumber.(eachrow(DP))

    DS = stack(DataFrame(DS, :auto))
    rename!(DS, "variable" => "time", "value" => "DS")
    select!(DS, Not(:time))
    DS.row = rownumber.(eachrow(DS))

    PP = stack(DataFrame(PP, :auto))
    rename!(PP, "variable" => "time", "value" => "PP")
    select!(PP, Not(:time))
    PP.row = rownumber.(eachrow(PP))

    SH = stack(DataFrame(SH, :auto))
    rename!(SH, "variable" => "time", "value" => "SH")
    select!(SH, Not(:time))
    SH.row = rownumber.(eachrow(SH))

    DF = stack(DataFrame(DF, :auto))
    rename!(DF, "variable" => "time", "value" => "DF")
    select!(DF, Not(:time))
    DF.row = rownumber.(eachrow(DF))

    df =  innerjoin(YT, CT, on = :row)
    df =  innerjoin(df, SH, on = :row)
    df =  innerjoin(df, HT, on = :row)
    df =  innerjoin(df, ST, on = :row)
    df =  innerjoin(df, LT, on = :row)
    df =  innerjoin(df, VT, on = :row)
    df =  innerjoin(df, PR, on = :row)
    df =  innerjoin(df, OT, on = :row)
    df =  innerjoin(df, DT, on = :row)
    df =  innerjoin(df, DS, on = :row)
    df =  innerjoin(df, DP, on = :row)
    df =  innerjoin(df, PP, on = :row)
    df =  innerjoin(df, DF, on = :row)
    df = innerjoin(df, ptmat, on = :time)

    if save == true
        CSV.write(raw"C:\Users\gdimeo\polybox\Dokumente\first_chapter\Model\results\mortgage_grid\tests\df_d.csv", df)
    end
    return(df)
end

function data_prep(df)
    # Generate an id by grouping by `time` and using row number within each group
    transform!(groupby(df, [:time]), :time => eachindex => :id)

    # Pull id to the fron 
    select!(df, [:id, :time], :)

    # Sort by id and time
    sort!(df, [:id, :time])

    df.time .= parse.(Int64, df.time)

    df.age .= df.time .+ 24
    # Sort by id and time
    sort!(df, [:id, :time])

    return(df)
end 
###########################################################################################################
begin
    C_E_s = zeros(Float64,(T, points, points, points, length(h_bar)))
    V_E_s = zeros(Float64,(T, points, points, points, length(h_bar)))

    C_E_r = zeros(Float64,(T, points, points))
    H_E_r = zeros(Float64,(T, points, points))
    V_E_r = zeros(Float64,(T, points, points))

    C_E_o = zeros(Float64,(T, points, points, points, length(h_bar)))
    V_E_o = zeros(Float64,(T, points, points, points, length(h_bar)))

    C_E_p = zeros(Float64,(T, points, points))
    H_E_p = zeros(Float64,(T, points, points))
    V_E_p = zeros(Float64,(T, points, points))
    print("")
end

Random.seed!(42); vect_sh = rand(2:N, 10000)
Random.seed!(42); YT_init = rand(LogNormal(lognormal_mu,lognormal_st), 10000) .* 12 ./ 10000 # draw  a starting incomes from lognormal distribution
Random.seed!(42); rn = (1 .+ rand(Normal(0, 0.1,), 10000*N))

function iter(mpb, val)
    θ_iter = (mpb*R/(1-mpb))^(γ)/(β*R)
    k_iter =  val*(β*R*θ_iter)^(1/γ)

    # 1: Older Renters
    begin  
        for t = 0:0
            #println(t)
            for y in 1:points
                result = dcegm_r(S[:], y, t, k_iter, θ_iter)
                C_E_r[T-t, y, :] = result[1]  

                H_E_r[T-t, y, :] = result[2] 

                V_E_r[T-t, y, :] = result[3]      
            end
        end
    end

    # 2: Older owners 
    begin #Owners
        for t = 0:0
            #println(t)
            for y in 1:points
                for o in 1:lastindex(ω)
                    for h in 1:lastindex(h_bar)
                        #println(y, "  ,  ",  o, " , ",  h, " , ", t)
                        result = dcegm_o(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_o[T-t, y, :, o, h] = result[1] 

                        V_E_o[T-t, y, :, o, h] = result[2]   
                    end
                end
            end
        end
    end

    # All other periods
    begin  
        for t = 1:T-1
            #println(t)
            for y in 1:points
                for o in 1:lastindex(ω)
                    for h in 1:lastindex(h_bar)
                        #println(y, " , ", " , ",  o, " , ",  h, " , ", t, k_iter, θ_iter)
                        result = dcegm_s(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_s[T-t, y, :, o, h] = result[1] 

                        V_E_s[T-t, y, :, o, h] = result[2] 
                        
                        result = dcegm_o(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_o[T-t, y, :, o, h] = result[1] 

                        V_E_o[T-t, y, :, o, h] = result[2]   
                    end 
                end
                result = dcegm_p(S[:], y, t, k_iter, θ_iter)
                C_E_p[T-t, y, :] = result[1]  

                H_E_p[T-t, y, :] = result[2] 

                V_E_p[T-t, y, :] = result[3]   

                result = dcegm_r(S[:], y, t, k_iter, θ_iter)
                C_E_r[T-t, y, :] = result[1]  

                H_E_r[T-t, y, :] = result[2] 

                V_E_r[T-t, y, :] = result[3]    
            end
        end
    end

    df = df_creation(10000, 10000, 0, 0, false)
    df_main = data_prep(df)
    df_main = select(@rsubset(df_main, :time <= 70), [:id, :time, :age, :DT])
    
    # Polynomial in simulated data
    lm1 = lm(@formula(DT ~ age + age^2  +  age^3 ), df_main)
    ages = collect(25:1:95)
    newdata_simul = DataFrame(age = ages)
    newdata_simul.age2 = newdata_simul.age .^ 2
    newdata_simul.age3 = newdata_simul.age .^ 3
    predicted_DT_simul = predict(lm1, newdata_simul, interval=:confidence)
    predicted_DT_simul.mpb .= mpb
    predicted_DT_simul.val .= val

    return(predicted_DT_simul)
end

mpb_val =  range(0.65,0.9,30)
val_val = range(0.1,3,30)

cd(results_smm)
for mv = 1:15 #lastindex(mpb_val)  
    for vv = 1:30 #lastindex(val_val)
            println("Processing mv=$mv, vv=$vv")
        try        
            df = iter(mpb_val[mv], val_val[vv])
            name = Dates.format(Dates.now(), "ddmmyy_HHMMSS")
            CSV.write(name * ".csv", df)

          catch e
            @warn "Error at mv=$mv, vv=$vv: $e"
            continue
        end
    end 
end

cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))

shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
shp = select(shp, [:id, :time, :age, :DT, :YT])
shp = @rsubset(shp, :DT > -1)
# Polynomial in SHP ()
lm1 = lm(@formula(DT ~ age + age^2  +  age^3 ), shp)
ages = collect(25:1:95)
newdata_shp = DataFrame(age = ages)
newdata_shp.age2 = newdata_shp.age .^ 2
newdata_shp.age3 = newdata_shp.age .^ 3
predicted_DT_shp = predict(lm1, newdata_shp, interval=:confidence)

cd(results_smm)
list_files = readdir(results_smm)
least_squares = DataFrame(values = Float64[], mpb = Float64[], val = Float64[])
for l = 1:lastindex(list_files)
    df = CSV.read(list_files[l], DataFrame)
    df.shp = predicted_DT_shp.prediction
    df.diff = (df.prediction .- df.shp) .^ 2
    push!(least_squares, [sum(df.diff)  df.mpb[1]  df.val[1]])
end

sort!(least_squares, [:mpb, :val])

z = fill(NaN, size(mpb_val, 1), size(val_val, 1))

for j in 1:lastindex(val_val)
    for i = 1:lastindex(mpb_val)
        z[i, j] = least_squares.values[i + lastindex(mpb_val) * (j -1)]
        #if  z[i, j] > 2.5
        #    z[i, j] = 2.5
        #end
    end
end

heatmap(mpb_val,
        val_val,
        z,
    xlabel="MPB", ylabel="10,000 CHF")


cd(results)
using Measures
plot!(size=(500,500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig(raw"smm.pdf")


#####
# Generate data with optimal calibration 
begin
    C_E_s = zeros(Float64,(T, points, points, points, length(h_bar)))
    V_E_s = zeros(Float64,(T, points, points, points, length(h_bar)))

    C_E_r = zeros(Float64,(T, points, points))
    H_E_r = zeros(Float64,(T, points, points))
    V_E_r = zeros(Float64,(T, points, points))

    C_E_o = zeros(Float64,(T, points, points, points, length(h_bar)))
    V_E_o = zeros(Float64,(T, points, points, points, length(h_bar)))

    C_E_p = zeros(Float64,(T, points, points))
    H_E_p = zeros(Float64,(T, points, points))
    V_E_p = zeros(Float64,(T, points, points))
    print("")
end


function iter2(mpb, val)
    θ_iter = (mpb*R/(1-mpb))^(γ)/(β*R)
    k_iter =  val*(β*R*θ_iter)^(1/γ)

    # 1: Older Renters
    begin  
        for t = 0:0
            #println(t)
            for y in 1:points
                result = dcegm_r(S[:], y, t, k_iter, θ_iter)
                C_E_r[T-t, y, :] = result[1]  

                H_E_r[T-t, y, :] = result[2] 

                V_E_r[T-t, y, :] = result[3]      
            end
        end
    end

    # 2: Older owners 
    begin #Owners
        for t = 0:0
            #println(t)
            for y in 1:points
                for o in 1:lastindex(ω)
                    for h in 1:lastindex(h_bar)
                        #println(y, "  ,  ",  o, " , ",  h, " , ", t)
                        result = dcegm_o(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_o[T-t, y, :, o, h] = result[1] 

                        V_E_o[T-t, y, :, o, h] = result[2]   
                    end
                end
            end
        end
    end

    # All other periods
    begin  
        for t = 1:T-1
            #println(t)
            for y in 1:points
                for o in 1:lastindex(ω)
                    for h in 1:lastindex(h_bar)
                        #println(y, " , ", " , ",  o, " , ",  h, " , ", t, k_iter, θ_iter)
                        result = dcegm_s(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_s[T-t, y, :, o, h] = result[1] 

                        V_E_s[T-t, y, :, o, h] = result[2] 
                        
                        result = dcegm_o(S[:], y, o, h, t, k_iter, θ_iter)
                        C_E_o[T-t, y, :, o, h] = result[1] 

                        V_E_o[T-t, y, :, o, h] = result[2]   
                    end 
                end
                result = dcegm_p(S[:], y, t, k_iter, θ_iter)
                C_E_p[T-t, y, :] = result[1]  

                H_E_p[T-t, y, :] = result[2] 

                V_E_p[T-t, y, :] = result[3]   

                result = dcegm_r(S[:], y, t, k_iter, θ_iter)
                C_E_r[T-t, y, :] = result[1]  

                H_E_r[T-t, y, :] = result[2] 

                V_E_r[T-t, y, :] = result[3]    
            end
        end
    end

    df = df_creation(10000, 10000, 0, 0, false)
    df_s = df_creation(10000, 10000, 1, 0, false)

    df_main = data_prep(df)
    df_perm = data_prep(df_s)

    return(df_main, df_perm)
end


function winsor2(var, prop)
    per = percentile.(eachcol(var), prop)
    var = ifelse.(var .> per, per, var)
    return(var)
end


least_squares[findmin(least_squares.values)[2], 2]
least_squares[findmin(least_squares.values)[2], 3]
val_theta = (least_squares[findmin(least_squares.values)[2], 2]*R/(1-least_squares[findmin(least_squares.values)[2], 2]))^(γ)/(β*R)
least_squares[findmin(least_squares.values)[2], 3]*(β*R*val_theta)^(1/γ)


dfs = iter2(least_squares[findmin(least_squares.values)[2], 2], least_squares[findmin(least_squares.values)[2], 3])

df_main = copy(dfs[1])
df_counterfactual1 = copy(dfs[2])

#####
default(;fontfamily="serif-roman")
cd(background_data)
index_df = CSV.read("index.csv", DataFrame)
index_df.time .= index_df.year .- 1977
select!(index_df, [:time, :year, :index_forecast])
leftjoin!(df_main, index_df, on = :time) 
df_main.tot_wealth = df_main.ST .+ df_main.DT .* df_main.HT .* df_main.pt .* df_main.OT
df_main.property = df_main.DT .* df_main.HT .* df_main.pt .* df_main.OT
df_main.ST .= df_main.ST ./ df_main.index_forecast ./ 10
df_main.tot_wealth .= df_main.tot_wealth ./ df_main.index_forecast ./ 10
df_main.YT .= df_main.YT ./ df_main.index_forecast ./ 10
df_main.property = df_main.property ./ df_main.index_forecast ./ 10
df_main.ST = winsor2(df_main.ST, 95)
df_main.tot_wealth = winsor2(df_main.tot_wealth, 95)
df_main.property = winsor2(df_main.property, 95)
df_main.data  .= "simul"
df_main.age .= df_main.time .+ 24
df_main = select(@rsubset(df_main, :time <= 70), [:id, :time, :age, :DT, :YT])
df_mean_simul = combine(groupby(df_main, [:age]), [:DT, :YT] .=> mean, renamecols = false)
df_mean_simul.data  .= "simul"

cd(shp_dir)
shp = DataFrame(XLSX.readtable("shp_allwaves_ready_hh.xlsx", "Sheet1"))
shp.DT .= shp.owner 
shp.id .= shp.idhous 
shp.age .= shp.avg_age
shp.time .= shp.year
shp.YT .=  ifelse.((shp.single .== 1), shp.income_hh, shp.income_hh./2)
shp.YT .= shp.YT ./ 100000
shp = select(shp, [:id, :time, :age, :DT, :YT])
shp = @rsubset(shp, :DT > -1)
shp = combine(groupby(shp, [:age]), [:DT, :YT] .=> mean, renamecols = false)
shp.data  .= "shp"

# Polynomial in SHP ()
lm1 = lm(@formula(YT ~ age + age^2  +  age^3 ), shp)
ages = collect(25:1:95)
newdata_shp = DataFrame(age = ages)
newdata_shp.age2 = newdata_shp.age .^ 2
newdata_shp.age3 = newdata_shp.age .^ 3

predicted_YT_shp = predict(lm1, newdata_shp, interval=:confidence)

# Polynomial in simulated data
lm1 = lm(@formula(YT ~ age + age^2  +  age^3 ), df_main)
newdata_simul = DataFrame(age = ages)
newdata_simul.age2 = newdata_simul.age .^ 2
newdata_simul.age3 = newdata_simul.age .^ 3

predicted_YT_simul = predict(lm1, newdata_simul, interval=:confidence)

plot(ages, 
predicted_YT_simul.prediction,
ribbon = (predicted_YT_simul.upper - predicted_YT_simul.prediction, predicted_YT_simul.prediction - predicted_YT_simul.lower),
title = "A: Per-capita household income",
xlabel = "Age",
ylabel = "100,000 CHF",
color = "#EA3680",
label="Simulation",
xticks = 25:10:95,
ylimits = (0,1),
legend = :topleft,
legend_font_pointsize = 11)

plot!(ages, 
predicted_YT_shp.prediction,
ribbon = (predicted_YT_shp.upper - predicted_YT_shp.prediction, predicted_YT_shp.prediction - predicted_YT_shp.lower),
color = "navy",
label="SHP",
linestyle = :dashdot)

p1 =plot!()

###
lm1 = lm(@formula(DT ~ age + age^2  +  age^3 ), shp)
ages = collect(25:1:95)
newdata_shp = DataFrame(age = ages)
newdata_shp.age2 = newdata_shp.age .^ 2
newdata_shp.age3 = newdata_shp.age .^ 3

predicted_DT_shp = predict(lm1, newdata_shp, interval=:confidence)

# Polynomial in SHP (ownership) 
lm1 = lm(@formula(DT ~ age + age^2  +  age^3 ), df_main)
ages = collect(25:1:95)
newdata_simul = DataFrame(age = ages)
newdata_simul.age2 = newdata_simul.age .^ 2
newdata_simul.age3 = newdata_simul.age .^ 3

predicted_DT_simul = predict(lm1, newdata_simul, interval=:confidence)

plot(ages, 
predicted_DT_simul.prediction,
ribbon = (predicted_DT_simul.upper - predicted_DT_simul.prediction, predicted_DT_simul.prediction - predicted_DT_simul.lower),
title = "B: Share of Homeowners",
xlabel = "Age",
ylabel = "Share",
color = "#EA3680",
label="Simulation",
xticks = 25:10:95,
ylimits = (0,1),
legend = :topleft,
legend_font_pointsize = 11)

plot!(ages, 
predicted_DT_shp.prediction,
ribbon = (predicted_DT_shp.upper - predicted_DT_shp.prediction, predicted_DT_shp.prediction - predicted_DT_shp.lower),
color = "navy",
label="SHP",
linestyle = :dashdot)

p2 =plot!()

cd(results)
plot(p1, p2, layour = (1,1))
using Measures
plot!(size=(1000,500), left_margin=7mm, 
    bottom_margin = 8mm, top_margin = 7mm) 
savefig(raw"shp_smoothed.pdf")



#####
# No transaction costs
sell_costs = 0
sell_costs1 = 0

dfs_app = iter2(least_squares[findmin(least_squares.values)[2], 2], least_squares[findmin(least_squares.values)[2], 3])

df_main_app = copy(dfs_app[1])
df_counterfactual1_app = copy(dfs_app[2])
