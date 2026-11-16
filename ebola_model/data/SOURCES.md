# Data sources

Both files in this folder were produced from the raw downloads kept in `data_raw/`. Downloaded on 16 November 2026.

## kikwit_1995_daily.csv

- **Content:** daily new cases by date of symptom onset and daily deaths, Kikwit, Democratic Republic of the Congo, 6 January to 16 July 1995. 192 rows, 292 onsets, 236 deaths.
- **Original source:** Khan AS, et al. The reemergence of Ebola hemorrhagic fever, Democratic Republic of the Congo, 1995. *J Infect Dis* 1999; 179 (Suppl 1): S76-S86.
- **Obtained from:** dataset `ebola_kikwit_1995` in the `outbreaks` R package (R Epidemics Consortium). Data provided to the package by T.J. McKinley, documented by Simon Frost.
  https://raw.githubusercontent.com/reconhub/outbreaks/master/data/ebola_kikwit_1995.RData
- **Licence of the package:** MIT.
- **Processing:** converted from RData to CSV with `tools/convert_kikwit.py`; columns renamed; a `day` column added (day 0 = 6 January 1995). No values changed.
- **Columns:** `date`, `day`, `onsets`, `deaths`, `daily_reporting` (False where the source did not record daily; all before 1 March 1995).
- **Known gaps:** 315 cases were identified in total; onset dates are missing for the rest. Recovery dates are not recorded.

## west_africa_2014_cumulative.csv

- **Content:** cumulative reported Ebola cases and deaths for Guinea, Liberia and Sierra Leone, 22 March 2014 to 5 January 2015. 263 rows.
- **Original source:** World Health Organization Disease Outbreak News and WHO Ebola situation reports.
- **Obtained from:** `country_timeseries.csv` in "Data for the 2014 Global Ebola outbreak", compiled by Caitlin Rivers.
  https://raw.githubusercontent.com/cmrivers/ebola/master/country_timeseries.csv
- **Licence of the compilation:** copying, redistribution, derivative works and commercial use are allowed; the source must be credited; it may not be sublicensed or redistributed under a different licence.
- **Processing:** reshaped from one column per country to one row per date and country; rows with no case count dropped. No values changed.
- **Columns:** `date`, `country`, `cumulative_cases`, `cumulative_deaths`.
- **Known problems:** counts are as reported at the time, were revised often (sometimes downwards) and under-count true cases. This is a third-party transcription and has not been checked line by line against the WHO documents.
