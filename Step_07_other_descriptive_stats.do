global path_simul "~results"

cd $path_simul
import delimited "dataready.csv", clear

keep id time ht ot pt yt lt dt ds dp event sh counterfactual 

gen flagdt = 1 if event == -1 & dt == 1
bys id: ereplace flagdt = max(flagdt)
keep if flagdt == 1
drop flagdt

keep if event >= 0

gen flagds = event if ds == 1 /* mark period of house sale */
bys id (time): ereplace  flagds = min(flagds) /* take first period */

gen censor = flagds != .

replace flagds = 30 if flagds == . 

bys id (time): gen dp1 = dp[_n+1]
sum dp1 if event == flagds & censor == 1 & counterfactual == 0 /* how many agents buy after selling */
sum dp1 if event == flagds & censor == 1 & counterfactual == 1

bys id (time): gen ht1 = ht[_n+1]
gen hdiff =  ht1 - ht 

sum hdiff if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1, d
sum hdiff if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1, d

ttest hdiff == 0 if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1
ttest hdiff == 0 if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1

gen hdiff_rel =  ht1/ht - 1

sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1
sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1



count if dp1 == 1  & event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1
count if dp1 == 1 &  event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1

count if hdiff < 0 & event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1
count if  hdiff < 0 &  event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1


* Effect heterogeneity - age
sum dp1 if event == flagds & censor == 1 & counterfactual == 0 & sh <= 20 /* how many agents buy after selling */
sum dp1 if event == flagds & censor == 1 & counterfactual == 1  & sh <= 20

sum hdiff if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & sh <= 20 , d
sum hdiff if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & sh <= 20 , d

sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & sh <= 20
sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & sh <= 20


sum dp1 if event == flagds & censor == 1 & counterfactual == 0 & sh > 20 /* how many agents buy after selling */
sum dp1 if event == flagds & censor == 1 & counterfactual == 1  & sh > 20

sum hdiff if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & sh > 20 , d
sum hdiff if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & sh > 20 , d

sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & sh > 20
sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & sh > 20


* Effect heterogeneity - income 
cd $path_simul
import delimited "dataready.csv", clear

qui sum yt if event == - 1, d
sca md = r(p50)

gen income_flag = 0 if event == -1 & yt <= md
replace income_flag = 1 if event == -1 & yt > md
bys id: ereplace income_flag = max(income_flag)

gen flagdt = 1 if event == -1 & dt == 1
bys id: ereplace flagdt = max(flagdt)
keep if flagdt == 1
drop flagdt

keep if event >= 0

gen flagds = event if ds == 1 /* mark period of house sale */
bys id (time): ereplace  flagds = min(flagds) /* take first period */

gen censor = flagds != .

replace flagds = 30 if flagds == . 

bys id (time): gen dp1 = dp[_n+1]
bys id (time): gen ht1 = ht[_n+1]
gen hdiff =  ht1 - ht 
gen hdiff_rel =  ht1/ht - 1


sum dp1 if event == flagds & censor == 1 & counterfactual == 0 & income_flag == 0 /* how many agents buy after selling */
sum dp1 if event == flagds & censor == 1 & counterfactual == 1  & income_flag == 0

sum hdiff if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & income_flag == 0
sum hdiff if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & income_flag == 0

sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & income_flag == 0
sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & income_flag == 0



sum dp1 if event == flagds & censor == 1 & counterfactual == 0 & income_flag == 1 /* how many agents buy after selling */
sum dp1 if event == flagds & censor == 1 & counterfactual == 1  & income_flag == 1

sum hdiff if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & income_flag == 1
sum hdiff if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & income_flag == 1

sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 0 & ds == 1 & dp1 == 1 & income_flag == 1
sum hdiff_rel if event == flagds & censor == 1 & counterfactual == 1 & ds == 1 & dp1 == 1 & income_flag == 1


sum ht if dt == 0 & counterfactual == 0 & event == 0 & income_flag == 1, d
sum ht if dt == 1 & counterfactual == 0 & event == 0 & income_flag == 1, d

sum h_costs if dt == 0 & counterfactual == 0 & event == 0 & income_flag == 1, d
sum h_costs if dt == 1 & counterfactual == 0 & event == 0 & income_flag == 1, d


***********************

cd $path_simul
import delimited "dataready.csv", clear


keep id time event sh dt ht yt counterfactual

replace id = id - 20000 if  counterfactual == 1

reshape wide dt ht yt, i(id time event) j(counterfactual)


preserve
collapse (mean) dt0 dt1 yt0 yt1, by(event)

gen change = ((dt1 - dt0)/dt0)/((yt1 - yt0)/yt0)

line change event

sum change if event >= 0 & event <= 30, d
restore 



preserve
keep if time > 20 & time <= 40 
collapse (mean) dt0 dt1 yt0 yt1, by(event)

gen change = ((dt1 - dt0)/dt0)/((yt1 - yt0)/yt0)

sum change if event >= 0 & event <= 30, d
restore 



preserve
replace ht0 = 0 if dt0 == 0 
replace ht1 = 0 if dt1 == 0 

collapse (mean) ht0 ht1 yt0 yt1, by(event)

gen change = ((ht1 - ht0)/ht0)/((yt1 - yt0)/yt0)

sum change if event >= 0 & event <= 30, d
restore 


***********************

cd $path_simul
import delimited "dataready.csv", clear

preserve
keep id event counterfactual dt 


collapse (mean) dt, by(event counterfactual)

sum dt if counterfactual == 0 & event == 30, d
sum dt if counterfactual == 1 & event == 30, d

restore 


cd $path_simul
import delimited "dataready_app.csv", clear

preserve
keep id event counterfactual dt 


collapse (mean) dt, by(event counterfactual)

sum dt if counterfactual == 0 & event == 30, d
sum dt if counterfactual == 1 & event == 30, d

restore 





