# exoplanet-characterization
Determines the radius and relative time of mid-transit of an exoplanet

The transit pipeline converts calibrated relative flux values from a given exoplanet's host star into a normalized and modeled transit light curve graph. These relative flux values are based on data that have been generated through the AstroImageJ graphical user interface.

The objective of this program is to create and model the normalized transit light curve graph of any exoplanet, given the light flux values and the corresponding Julian Dates of each calibrated host-star field image. The modeled transit light curve graphs are then utilized to calculate the planet's normalized radii (Rp/R*) and relative time of mid-transit. Chi-squared goodness of fit tests are performed on the data to find the best-fit model. Additionally, chi-squared maps are used to find the 1σ error bar ranges for all of the light curve calculations.

Created by Sarah Tang on Wed October 24 16:03:10 2018 with assistance from William Waalkes.
