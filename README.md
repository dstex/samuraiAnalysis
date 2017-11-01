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
* [simplekml](https://simplekml.readthedocs.io/en/latest/) : _Used only for generating KML/KMZ files for use in Google Earth_

## Functionality
samuraiAnalysis is largely just a collection of functions designed to simplify the plotting of SAMURAI output. The following is a list of all modules and included functions, complete with brief descriptions of capabilities and use.
