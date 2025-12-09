[README_MAIN_ChANGE_PulsePress_SGS.txt](https://github.com/user-attachments/files/24063243/README_MAIN_ChANGE_PulsePress_SGS.txt)
####
####
####
####Data from: Assessing interactions between deluge pulse and nitrogen press on productivity in a semi-arid grassland: the importance of timing

[Access this dataset on Dryad](DOI link in process)



Persistent increases in nitrogen availability represent a chronic resource “press.” Concurrently, intensifying weather can create more extreme rainfall events, or "deluges," which represent a resource “pulse.” Previous experiments that assess short-term alteration of nitrogen and water often find additive effects, but theory and the pulse-press framework suggest that a single extreme deluge amid the backdrop of chronic nitrogen addition may result in interaction. Synergistic interaction between nitrogen press and deluge pulse would result in combined effects greater than the sum of their separate parts, but explicit tests of pulse-press dynamics are exceedingly rare. Given that the magnitude of nitrogen differentially affects community structure and function, the amount of nitrogen addition may also mediate interaction between nitrogen and deluge. In this experiment, we leveraged a long-term nitrogen addition study and two years of experimental deluge to create two case studies that enable us to assess the role of role of nitrogen magnitude in moderating pulse-press interactions. 

This data originates from a shortgrass steppe site of northeastern Colorado in which nitrogen (as urea) has been added since 2014 at eight levels: 0, 2.5, 5, 7.5, 10, 15, 20, 30 g m-2. This experiment used a randomized, complete block design with six blocks composed of eight plots each receiving one of the eight nitrogen addition levels. Deluge subplots were nested within each plot. Experimental deluge events represented the 95th percentile of historic precipitation events. One deluge event was applied mid-late season (late July to early August) in 2021 (during the eighth year of nitrogen addition) and a second was applied early-mid season (late June) in 2022 (during the ninth year of nitrogen addition) to separate subplots within each nitrogen addition plot.

The 2021 mid-late season deluge had no significant effect on ANPP and did not interact with nitrogen addition. In contrast, the 2022 early-mid deluge increased ANPP and synergistically interacted with nitrogen addition. The ANPP response with deluge was unimodal across the nitrogen gradient, with an estimated peak at ~17 g m-2 of nitrogen. The 2022 synergistic interaction appeared driven by fast-growing, weedy forbs, which increased in abundance from 2021 to 2022 and were positively effected by nitrogen addition. Forb biomass in 2022 increased by 16% across nitrogen addition levels, 206% in response to deluge, and 1088% in response to nitrogen and deluge. Species richness and evenness were effected by nitrogen and deluge, but not their interaction. These results demonstrate that long-term nitrogen addition and short-term deluge can synergistically interact to effect a key ecosystem response and that the magnitude of nitrogen addition can mediate the magnitude of synergy, however these effects may be strongly moderated by environmental context.

-----Full details and methodology can be found in the associated manuscript-----



###
###
###Data description --

There are six files associated with this project, each with a single major type of data: 1. ANPP, 2. canopy greenness, 3. ammonium and nitrate availability, 4. precipitation, 5. soil moisture, 6. plant community composition.

###1.
#MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022:
	ANPP was collected at the end of each growing season at peak growth in 2021 and 2022 from all nitrogen plots and deluge subplots. The ANPP sampling area was never resampled and separate from the permanent species composition plots. ANPP was estimated by clipping plant biomass down to the ground level within two 20 x 50 cm quadrats. Clipped plants were sorted into functional groups: grass, forb, woody. Biomass was then dried, dead biomass was sorted out, and the samples were weighed. Final ANPP is calculated as the sum of forb and grass biomass, multiplied by 10 to achieve ANPP g m-2. Woody plants were uncommon and not representative of the site, so woody biomass was excluded from the final ANPP measurement.

#DATAFRAME COLUMNS:
block: 6 blocks from A-F
plot: 48 plots from 1-48
site: site name "SGS," an abbreviation for "shortgrass steppe"
year: year of data collection, 2021 or 2022
treatment: three levels--N (nitrogen-only), PP1 (pulse-press plot where deluge was applied in 2021), PP2 (pulse-press plot where deluge was applied in 2022)
grass: g m-2 of grass biomass
forb: g m-2 of forb biomass
woody: g m-2 of woody biomass
dead: g m-2 of dead biomass that was sorted out of the sample prior to final weighing
cactus: The proportion of cactus within the each sampling area
anpp: g m-2 of grass and forb biomass
anpp_with_woody: g m-2 of grass, forb, and woody biomass
nitrogen: g m-2 of nitrogen added to plots, eight levels--0, 2.5, 5, 7.5, 10, 15, 20, 30


###2.
#MAIN_ChANGE_PulsePress_SGS_calculated green up values_2021to2022:
	Canopy greenness data was collected before deluge application until the end of the growing season in all nitrogen plots and deluge subplots through repeat photography. The green chromic coordinate (GCC) of each photo was calculated as the ratio of green to red and blue light within each pixel then averaged across all pixels in each image.

#DATAFRAME COLUMNS:
block: 6 blocks from a-f
site: site name "sgs," an abbreviation for "shortgrass steppe"
year: year of data collection, 2021 or 2022
plot: 48 plots from 1-48
treatment: three levels--n (nitrogen-only), pp1 (pulse-press plot where deluge was applied in 2021), pp2 (pulse-press plot where deluge was applied in 2022)
date: the date data was collected formatted as mm/dd/yyyy
gcc_mean: the mean GCC 
doy: day of year when data was collected
dom: day of month data was collected
month: the month of data collection
nitrogen: g m-2 of nitrogen added to plots, eight levels--0, 2.5, 5, 7.5, 10, 15, 20, 30
deluge.date.1: the day the first-half of experimental deluge was applied, formatted as mm/dd/yyyy
deluge.date.2: the day the second-half of experimental deluge was applied, formatted as mm/dd/yyyy


###3.
#MAIN_ChANGE_PulsePress_SGS_corrected ammonium nitrate_2021to2022:
	Biological available nitrate and ammonium data were collected using mixed bed ion exchange resin bags buried in nitrogen plots and corresponding deluge subplots. Resin bags were buried during the period of experimental deluge and removed when soil moisture no longer differed between nitrogen plots and deluge subplots. Resin bags were extracted with KCL then samples were analyzed on a O.I. Analytical 3700 Automated Chemistry Analyzer to assess ammonium and nitrate content.

#DATAFRAME COLUMNS:
site: site name "sgs," an abbreviation for "shortgrass steppe"
sample_year: year of data collection, 2021 or 2022
plot: block (A-F) and plot (1-48)
nitrogen: g m-2 of nitrogen added to plots, eight levels--0, 2.5, 5, 7.5, 10, 15, 20, 30
treatment: three levels--n (nitrogen-only), pp1 (pulse-press plot where deluge was applied in 2021), pp2 (pulse-press plot where deluge was applied in 2022)
avgcalc_mg_l: the amount of ammonium or nitrate as mg l-1
channel: ammonium or nitrate
outlier: "x" indicates an outlier sample, or a sample associated with an outlier


###4.
#MAIN_ChANGE_PulsePress_SGS_precip_Jan1toAug31_2021to2022:
	Site-level daily precipitation was obtained from a nearby USDA precipitation gauge and can be found at: https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=2017. 

#DATAFRAME COLUMNS:
site: site name "sgs," an abbreviation for "shortgrass steppe"
year: year of data collection, 2021 or 2022
date: date of rainfall formatted as mm/dd/yyyy
precipitation_mm: the amount of incremental precipitation for each day in mm
precipitation_summed_mm: the total amount of precipitation received each day summed throughout the growing season in mm
doy: day of year when data was collected
dom: day of month data was collected
month: the month of data collection


###5.
#MAIN_ChANGE_PulsePress_SGS_soil moisture_2021to2022:
	Soil moisture,  measured as volumetric water content (% VWC), was collected before deluges then monitored until the end of the growing season using a HydroSense II Handheld Soil Moisture Sensor (Campbell Scientific, Inc., Logan, UT), equipped with 20 cm soil moisture probes. Measurements were taken in three random locations then averaged.

#DATAFRAME COLUMNS:
block: 6 blocks from a-f
site: site name "sgs," an abbreviation for "shortgrass steppe"
year: year of data collection, 2021 or 2022
plot: 48 plots from 1-48
treatment: three levels--nitrogen (nitrogen-only), pulse-press1 (pulse-press plot where deluge was applied in 2021), pulse-press2 (pulse-press plot where deluge was applied in 2022)
date: date of soil moisture collection  formatted as mm/dd/yyyy
vwc: volumetric water content (%) of soil
deluge.date.1: the day the first-half of experimental deluge was applied, formatted as mm/dd/yyyy
deluge.date.2: the day the second-half of experimental deluge was applied, formatted as mm/dd/yyyy
doy: day of year when data was collected
dom: day of month data was collected
month: the month of data collection
days_since_deluge: the amount of days since the final date of deluge application


###6.
#MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022:
	Species composition data was collected in early and late growing season to capture maximum coverage of plant species present. Abundance is collected as percent cover, visually estimated within a permanent 1 x 1 m quadrat.

#DATAFRAME COLUMNS:
block: 6 blocks from A-F
plot: 48 plots from 1-48
site: site name "sgs," an abbreviation for "shortgrass steppe"
treatment: three levels--nitrogen (nitrogen-only), pulse press 1 (pulse-press plot where deluge was applied in 2021), pulse press 2 (pulse-press plot where deluge was applied in 2022)
year: year of data collection, 2021 or 2022
species: species of plant, formatted as genus and specific epithet
cover: abundance of species as percent cover, integer value
nitrogen: g m-2 of nitrogen added to plots, eight levels--0, 2.5, 5, 7.5, 10, 15, 20, 30



###
###
###Sharing/Access information

Data and code can also be found in the personal GitHub account of the lead author: https://github.com/mlinabury/ChANGE-PulsePress-Manuscript-Files

Precipitation data was derived from a USDA rain gauge found at: https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=2017



###
###
###Code/Software

All analyses and data preparation were conducted within R, version 4.2.2.
