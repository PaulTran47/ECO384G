/******************************************************************************/
/* 1 asterisk represent commented out code */
/** 2 asterisks represent comments for code or commented out code **/

/* One asterisk represents commented out code. */
/** Two asterisks represent comments for code. **/

/* Answers to question parts that don't involve code can be found at the */
/* bottom of the programme, in the section ``Questions asked in problemset */
/* that don't involve code. */

/* Text answers to question parts that involve code will be between the */
/* sub-section label: */
/**********/
/* ANSWER */
/**********/
/* Answer here */
/**************/
/* END ANSWER */
/**************/

/* Comments that are important will be between the sub-section label: */
/********/
/* NOTE */
/********/
/* Important note here */
/************/
/* END NOTE */
/************/
/* ECO384G Problem Set 4 (1, Spring 2022;12, Nitya), 2 */
/* Paul Le Tran, plt377 */ 
/* 10 May, 2022 */
/******************************************************************************/

/******************************************************************************/
/** Setting up workspace **/
clear
cls

set matsize 547

local home_dir = "\\path\to\programmes"
local data_dir = "\\path\to\data"

cd `home_dir'
/******************************************************************************/

/******************************************************************************/
/** Part 2.1-2: Loading in data and setting up variables **/
cd `data_dir'
use gravity2000

/* Taking log of bilateral exports */
gen ln_biexports = log(biexports)

/* Taking log of distance */
gen ln_dist = log(dist)

/* Taking log of GDP variables */
gen ln_gdpnomusd1 = log(gdpnomusd1)
gen ln_gdpnomusd2 = log(gdpnomusd2)
gen ln_gdp12 = log(gdpnomusd1*gdpnomusd2)

cd `home_dir'
/******************************************************************************/

/******************************************************************************/
/** Part 2.3-5: Estimate the "traditional" gravity equation **/
regress ln_biexports ln_gdp12 ln_dist, robust cluster(dist)

/**********/
/* ANSWER */
/**********/
/* Given the results from our estimation of the traditional gravity equation, */
/* we see that if all else is held constant, distance between countries */
/* increasing by 1% is associated with bilateral exports decreasing by about */
/* 1.0171% on average. Therefore, in the situation that distance is doubled */
/* (i.e., increases by 100%), bilateral exports will decrease by -101.71% on */
/* average. */
/**************/
/* END ANSWER */
/**************/

/**********/
/* ANSWER */
/**********/
/* Given the results from our estimation of the traditional gravity equation, */
/* we see that if all else is held constant, the product of the countries' */
/* GDPs increasing by 1% is associated with bilateral exports increasing by */
/* roughly 0.8921% on average. Therefore, in the situation that said product */
/* doubles, bilateral exports will increase by 89.21% on average. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2.6: Adding extra variables to gravity equation **/
regress ln_biexports ln_gdp12 ln_dist i.border i.landl i.comlang_off i.colony, robust cluster(dist)

/**********/
/* ANSWER */
/**********/
/* Holding all else constant, we see that the effect of an increase by 1% in */
/* each of the following variables will be associated with the following per */
/* cent change in bilateral exports: Border = 0.8771%; */
/* Landlocked (1 country) = -0.4853%; Landlocked (2 countries) = -0.0650%; */
/* Common language = 0.6723%; Colony = 1.0382%. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2.7: If the underlying model of trade is Armington, is the traditional gravity equation misspecified, and why? **/
/**********/
/* ANSWER */
/**********/
/* The traditional gravity equation would indeed by misspecified because it */
/* would not be accounting for the effect of prices, which are relevant under */
/* the Armington trade model. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2.8: Include exporter and importer effects in the gravity equation. Can we now say anything about how country size (GDP) affects bilateral trade? **/
/********/
/* NOTE */
/********/
/* Not splitting by factors because we just care if country is an importer or */
/* an exporter, not if it has a unique code. */
/************/
/* END NOTE */
/************/
regress ln_biexports ln_gdp12 ln_dist i.ifscode1 i.ifscode2, robust cluster(dist)

/**********/
/* ANSWER */
/**********/
/* By including importer- and exporter-fixed effects, we see that the effect */
/* of ln(X_{i}X_{j}) on bilateral exports has increased from 0.8921% to */
/* 1.0061%. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2.9: How does doubling of distance change bilateral exports? Compare your answer to the traditional gravity estimate. **/
/**********/
/* ANSWER */
/**********/
/* By including importer- and exporter-fixed effects, we see that the effect */
/* of increasing distance by 1% on bilateral exports has changed from causing */
/* bilateral exports to decrease by -1.017% to -1.3550%. Therefore, doubling */
/* distance results in bilateral exports to decrease by 135.5% compared to */
/* old 101.7%. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2.10: Modified gravity equation with importer- and exporter-fixed effects **/
/********/
/* NOTE */
/********/
/* Not splitting by factors because we just care if country is an importer or */
/* an exporter, not if it has a unique code. */
/************/
/* END NOTE */
/************/
regress ln_biexports ln_gdp12 ln_dist i.border i.landl i.comlang_off i.colony i.ifscode1 i.ifscode2, robust cluster(dist)

/**********/
/* ANSWER */
/**********/
/* By including importer- and exporter-fixed effects, we see the following */
/* changes to the impact of increasing each of the following variables by 1% */
/* on bilateral exports in percentage terms: Border = 0.4167% from 0.8602%; */
/* Landlocked (1 country) = -1.4296% from -0.4853%; */
/* Landlocked (2 countries) = -2.7519% from -0.0650%; */
/* Common language = 0.5839 from 0.6723%; Colony = 0.9896% from 1.0382%. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/