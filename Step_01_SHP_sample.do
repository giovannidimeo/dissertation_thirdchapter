global path_shp "~" // path to folder wher to store data
global path_raw "~" // path to folder with raw SHP data
gloval path_background "~\background_data" // path to folder with background data


forval i = 0/20 { 
	cd $path_raw
	local k = `i' + 2
	if 	`i' < 10 {
		use "SHP-Data-W1-W24-STATA\W`k'_200`i'\shp0`i'_p_user.dta", clear // open individual file
		merge m:1 idhous0`i' using "SHP-Data-W1-W24-STATA\W`k'_200`i'\shp0`i'_h_user", nogen // add all variables from household file

		keep  idpers idhous0`i' age0`i' sex0`i' nbadul* wstat0`i' i0`i'ptotn h0`i'h29 h0`i'h06 h0`i'h07 cohast0`i' h0`i'h20  nbkid0`i' permit0`i' p0`i'i02 h0`i'i30
		keep if age0`i' >= 25 & age0`i' <= 95

		cap rename idhous0`i' idhous
		cap rename sex0`i' sex
		rename age0`i' age 
		cap rename nbadul0`i' nr_adults
		cap rename h0`i'h29 owner
		cap rename h0`i'h06 moved
		cap rename h0`i'h07 reason_moved /* 21, 21, 22 */
		cap rename wstat0`i' occup_stat
		cap rename i0`i'ptotn iptotn
		cap rename cohast0`i' cohast
		cap rename h0`i'h20  nr_rooms
		cap rename nbkid0`i'  nr_children
		cap rename permit0`i' permit
		cap rename p0`i'i02 change_finan
		cap rename h0`i'i30 satisf_finan

		gen wave = `i'

		cd $path_shp
		save "wave`i'", replace
	} 
	else {
		use "SHP-Data-W1-W24-STATA\W`k'_20`i'\shp`i'_p_user.dta", clear // open individual file
		merge m:1 idhous`i' using "SHP-Data-W1-W24-STATA\W`k'_20`i'\shp`i'_h_user", nogen // add all variables from household file

		keep  idpers idhous`i' age`i' sex`i' nbadul* wstat`i' i`i'ptotn h`i'h29 h`i'h06 h`i'h07 cohast`i' h`i'h20  nbkid`i' permit`i' p`i'i02 h`i'i30
		keep if age`i' >= 25 & age`i' <= 95

		cap rename idhous`i' idhous
		rename age`i' age 
		cap rename sex`i' sex
		cap rename nbadul`i' nr_adults
		cap rename h`i'h29 owner
		cap rename h`i'h06 moved
		cap rename h`i'h07 reason_moved /* 21, 21, 22 */
		cap rename wstat`i' occup_stat
		cap rename i`i'ptotn income_pers
		cap rename cohast`i' cohast
		cap rename h`i'h20  nr_rooms
		cap rename nbkid`i'  nr_children
		cap rename permit`i' permit
		cap rename p`i'i02 change_finan
		cap rename h`i'i30 satisf_finan

		gen wave = `i'

		cd $path_shp
		save "wave`i'", replace
	} 
} 


cd $path_shp
use wave0, clear
forval i = 1/20 { 
	append using wave`i'
}

gen year = 2000 + wave + 2
cd \\nash\mtec-home\gdimeo\Desktop\SHP\data\Data_STATA\Data_STATA
merge m:1 idhous year using "SHP-Data-Imputed-Income-Wealth-STATA\imputed_income_hh_long_shp", keep (1 3) nogen keepusing (ihtyni)
rename ihtyni income_hh

drop iptotn 
merge 1:1 idpers year using "SHP-Data-Imputed-Income-Wealth-STATA\imputed_income_pers_long_shp", keep (1 3) nogen keepusing (iptotn)


order idpers wave year 

sort idpers wave

save shp_allwaves, replace

global  $path_background
import delimited "C:\Users\gdimeo\polybox\Dokumente\third_chapter\index.csv", clear

gen time = year - 1977
drop index

gloval $path_shp
merge 1:m year using shp_allwaves, keep(3) nogen

foreach var in income_pers income_hh iptotni {
	replace `var' = `var'/index_forecast
}


gen owner2 = . 
replace owner2 = . if owner < 0 
replace owner2 = 0 if owner == 1 | owner == 3
replace owner2 = 1 if owner == 2
drop owner 
rename owner2 owner

cd $path_shp
save shp_allwaves, replace


******************************
cd $path_shp
use shp_allwaves, clear

keep if income_hh > 0 
keep if income_hh != . 
winsor2 income_hh iptotni, replace cuts(1 99)

* Restrict dataframe
keep idhous year sex age owner income_hh iptotni idpers nr_adults

* keep only 1- and 2-person households
sort year idhous
duplicates tag  idhous year, gen(flag)
keep if flag <= 1
drop flag 

* Keep only households younger than 60
bys idhous year: egen max_age = max(age)

* Reshape dataframe
drop idpers
sort idhous year sex
egen id_new = tag(idhous year) 

reshape wide iptotni sex age, i(idhous year) j(id_new)

gen single = 0 
replace single = 1 if sex0 == . | sex1 == . 

* Generate income change
bys idhous (year): gen hh_yt1 = income_hh[_n-1]
gen change_income = (income_hh)/(hh_yt1) - 1
keep if change_income != . 
winsor2 change_income, cuts(5 95) replace

* Standardize variables 
bys year: egen change_income_std = std(change_income)

keep if owner != . 

* Generate shock variable
cap drop flag* 

gen flag = change_income <= -0.3 
unique idhous if flag == 1

gen shock = year if flag == 1
bys idhous (year): ereplace shock = min(shock)
replace shock = 0 if mi(shock)
egen avg_age = rowmean(age0 age1)

gen event = year - shock if shock != 0 
sum change_income if shock != 0, d
sca std = r(sd)
gen flag2 =  abs(change_income) >= std/2
gen flag3 = .
replace flag3 = 1 if event == -1 & flag2 == 1
bys idhous: ereplace flag3 = max(flag3) 
rename flag3 increase_in_neg1
replace increase_in_neg1 = 0 if increase_in_neg1 == .
unique idhous if shock!= 0 & increase_in_neg1 == 0
drop flag2

cd $path_shp
save "shp_allwaves_ready_hh", replace

cd $path_shp
export excel using "shp_allwaves_ready_hh.xlsx", firstrow(variables) replace






