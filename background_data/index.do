import excel "C:\Users\gdimeo\polybox\Dokumente\first_chapter\Model\GitHub\background_data\consumer_prices.xlsx", sheet("Sheet1") clear
drop A	 
collapse (mean) C, by(B)
drop if C == .	
rename B year
rename C index
destring year, replace

set obs 95

forvalues i=1(1)48 {
	replace year[`k'] = year[47] + `i' if _n == 47 + `i'
}


ipolate index year, generate(index_forecast) epolate

* Save data
cd $main
save index.dta, replace
export delimited using "C:\Users\gdimeo\polybox\Dokumente\first_chapter\Model\GitHub\background_data\index.csv", replace
