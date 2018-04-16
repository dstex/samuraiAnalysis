# samuraiAnalysis
The samuraiAnalysis package is designed to simplify post-processing, plotting, and analysis of output from the Spline Analysis at Mesoscale Utilizing Radar and Aircraft Instrumentation ([SAMURAI](https://github.com/mmbell/samurai)) package. Currently, focus is placed on plotting, with post-processing and analysis methods to be added in the future.

## Dependencies
This package has been tested with both Python 2.7 and 3.5+. Installation and use of the [Anaconda](https://www.anaconda.com/distribution/) package is likely the easiest way to obtain most dependencies.
* [NumPy](http://www.scipy.org)
* [SciPy](http://www.scipy.org>)
* [matplotlib](http://matplotlib.org/)
* [Py-ART](https://github.com/ARM-DOE/pyart)
* [Cartopy](http://scitools.org.uk/cartopy/)
* [xarray](http://xarray.pydata.org/en/stable/)
* [pandas](http://pandas.pydata.org/)
* [FLplot](https://github.com/stechma2/FLplot/)
* [YAML](https://anaconda.org/anaconda/yaml/)
* [simplekml](https://simplekml.readthedocs.io/en/latest/) : _Used only for generating KML/KMZ files for use in Google Earth_

## IMPORTANT NOTE
### Cross-section location(s)
If specifying multiple cross-sections, format the `xsStrt`/`xsEnd` vars as lists of tuples in config file. For example:
```
xsStrt:
- !!python/tuple
  - 43.321
  - -97.841
- !!python/tuple
  - 43.2
  - -98.1
xsEnd:
- !!python/tuple
  - 44.214
  - -98.298
- !!python/tuple
  - 44.214
  - -98.298
```
Here, it is important to double check that the winds within the plane of the cross-section(s) match 
what is expected based on those seen in the PPI(s). If winds appear to be opposite of those expected, 
swap the initial and final XS points. _(Still trying to find a way to implement automatic checks for this.)_

## Functionality
samuraiAnalysis is largely just a collection of functions designed to simplify the plotting of SAMURAI output. A list of all modules and included 
functions, complete with brief descriptions of capabilities and use will be included here in the future.