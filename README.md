# MGM
Microdosimetric Gamma Model to calculate DNA damage

## Introduction
MGM is a microdosimetric gamma model to calculate DNA damage depending on the values of microdosimetric quantities for a given beam.

## Installation
### Requirements
* Python 3
* Python packages: numpy, scipy

## Usage
### Library
The library can be used to calculate the DNA damage provided yF for a monoenergetic beam:
```
import mgm
yF = 4.0 # microdosimetric quantity for a given beam in kev/um

# DNA damage for a monoenergetic beam
MGM = mgm.MicrodosimetricGammaModel()
base_damage = MGM.getBD(yF)
strand_breaks = MGM.getSB(yF)
direct_base_damage = MGM.getBDD(yF)
indirect_base_damage = MGM.getBDI(yF)
direct_strand_breaks = MGM.getSBD(yF)
indirect_strand_breaks = MGM.getSDI(yF)
num_sites = MGM.getN_sites(yF)
num_sites_with_DSB = MGM.getN_sites_with_DSB(yF)
complexity_distribution = MGM.getComplexityDistribution(yF)
```