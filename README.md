# exoplanets
Determines the radius and relative time of mid-transit of a given exoplanet.

The Exoplanet Transit Pipeline converts calibrated relative flux values from a given exoplanet's host star into a normalized and modeled transit light curve graph. These relative flux values are based on data that have been generated through the AstroImageJ graphical user interface.

The objective of this program is to create and model the normalized transit light curve graph of any exoplanet, given the light flux values and the corresponding Julian Dates of each calibrated host-star field image. The modeled transit light curve graphs are then utilized to calculate the planet's normalized radii (Rp/R*) and relative time of mid-transit. Chi-squared goodness of fit tests are performed on the data to find the best-fit model. Additionally, chi-squared maps are used to find the 1σ error bar ranges for all of the light curve calculations.

The exoplanet LHS 3844 b is analysed in the Exoplanet Transit Pipeline file in this repository. LHS 3844 b was an exoplanet detected by NASA's Transiting Exoplanet Survey Satellite in early 2020. I obtained raw data of the flux values of LHS 3844 b's host star during the time of the planet's transit from the El Sauce Observatory in Río Hurtado, Coquimbo, Chile.

Created by Sarah Tang on Wed October 24 16:03:10 2018 with assistance from William Waalkes.
