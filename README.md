# fetch-params-simbad
Scrip to fetch various stellar parameters from Simbad.


The main script 

> make-star-dataframe.py

is what you can run to get the parameters.
This script takes as an input, a text file containing a list of stellar designations, it matters not what catalogue or convention, since the script utilizes Simbad, and Simbad is intelligent. _**starlist.txt**_ is an example of such an input file.

The output of the main script is a csv file; _**star_data.csv**_.
