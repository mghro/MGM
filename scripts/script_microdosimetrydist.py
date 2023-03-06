#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/5/23 11:34 AM

@author: alejandrobertolet
"""

# Example of use
from src import mgm

# Data file, format is 'microdose' (3 columns with energy, specific energy and lineal energy)
data_file = 'xray_microdosimetry_1um.phsp'

# Other supported formats are:
# - A list of lineal energy values; each value in a new row
# - A list of bins with number of counts in each bin; each bin in a new row, separated by a space

# Initialize calculator.

calc = mgm.MicrodosimetryGammaCalculator(data_file, format='microdose', subsample=1000)
calc.CalculateDamage()
calc.PlotComplexityDistribution(density=True)

# Gets the number of sites with DSBs
print(calc.getNumberOfSitesWithDSB(perTrack=True))

# Distribute damages over a nucleus with radius 3.5 um. Dose scales the number of damage sites.
# This returns a list of damage sites with position (x, y, z), and complexity
dose = 2 #Gy
damages = calc.DistributeDamageOverNucleus(dose=dose, radius=3.5, inTracks=True)
print('Number of sites for ', dose, ' Gy: ', len(damages))
calc.PlotDistributedDamageOverNucleus()