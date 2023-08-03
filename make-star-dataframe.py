# -*- coding: utf-8 -*-
"""
Created: 2023-08-02
Last Edit: 2023-08-02
Author: Caleb Brand
Contact: brand.c.research@gmail.com
"""

# ############# IMPORTS SETUP

import argparse
import numpy as np
import pandas as pd
import datetime as dt
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Latitude

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, help="Relative path to the input file")
args = parser.parse_args()

inputfilepath = args.file or input("File path: ")


# ############# READ INPUT FILE
inputfiledata = pd.read_csv(
    inputfilepath,
    header=None,
    index_col=False,
    delimiter="\s+",
)


# ############# SET UP OUTPUT DATAFRAME
outputfiledata = pd.DataFrame(
    columns=["name",
        "hd",
        "entered_name",
        "ra",
        "dec",
        "type",
        "vmag",
        "bmv",
        "extinction",
        "teff",
        "logg",
        "feh",
        "visibility",
        "max_alt",
        "exp_time",
    ]
)


# ############# FUNKTION TO RETRIEVE DATA FROM ONLINE SOURCES

def estimate_exposure(magnitude):
    exposure_dictionary = {
        0 : 100,
        1 : 150,
        2 : 200,
        3 : 300,
        4 : 600,
        5 : 900,
        6 : 1200,
        7 : 1800,
    }
    rounded = np.round(magnitude)
    if rounded < 0:
        return 100
    elif rounded > 7:
        return 1800
    exposure = exposure_dictionary[rounded]
    return exposure



def retrieve_hd_number(designation):
    hd_query = Simbad.query_objectids(designation)
    # print(hd_query.items())
    for thing in hd_query:
        if "HD" in thing["ID"]:
            return thing["ID"]
    return designation



def generate_datetimes():
    start_date = dt.datetime(2023, 1, 1, 13, 0)
    end_date = dt.datetime(2023, 12, 31, 13, 0)
    # using 13:00 in the above two lines since 
    # we haven't corrected for LC = UT + 12
    datetimes = []
    current_date = start_date
    while current_date <= end_date:
        datetimes.append(current_date)
        current_date += dt.timedelta(days=1)
    return datetimes



def calculate_alt(designation, datetimes):
    radec = SkyCoord.from_name(designation)
    MJO = EarthLocation(
        lat=-43.9853*u.deg,
        lon=170.4641*u.deg,
        height=1029*u.m,
    )
    altazs = []
    for time in datetimes:
        altaz = radec.transform_to(AltAz(obstime=time, location=MJO))
        altazs.append(altaz.alt)
    max_alt, min_alt, alt_index = (max(altazs), min(altazs), altazs.index(max(altazs)))
    max_alt_date = datetimes[alt_index].date().month

    # we specifically choose 28 degrees in the following since that is 
    # around the lowest the McLellan telescope can observe in the south
    min_observable = Latitude(20, unit='deg')
    if min_alt >= min_observable:
        months_visible = "thoughout"
        return max_alt, max_alt_date, months_visible

    # VERY JANKY THING BELOW, PLEASE REWRITE BETTER
    months_visible = []
    for altitude in altazs:
        if altitude > min_observable:
            index = altazs.index(altitude)
            month = datetimes[index].date().month
            if month not in months_visible:
                months_visible.append(month)
    # YUCK YUCK YUCK

    return max_alt, max_alt_date, months_visible


def retrieve_data(designation, datetimes):
    # add fields to the default query
    Simbad.add_votable_fields(
        'flux(V)',
        'flux(B)',
        'fe_h',
        'sp',
    )
    # run the query
    main_query = Simbad.query_object(designation)
    # assign variables to some of the data for ease of processing
    hd_name = retrieve_hd_number(designation)
    vmag = main_query['FLUX_V'].item()
    exposure = estimate_exposure(vmag)
    max_alt, max_alt_date, months_visible = calculate_alt(designation, datetimes)
    # assign the query fields to a series for appending to the output file
    full_data = pd.Series(
        data = {
            "name" : main_query["MAIN_ID"].item(),
            "hd" : hd_name,
            "entered_name" : f"{designation}",
            "ra" : main_query["RA"].item(),
            "dec" : main_query["DEC"].item(),
            "type" : main_query["SP_TYPE"].item(),
            "vmag" : vmag,
            "bmv" : main_query['FLUX_B'].item() - main_query['FLUX_V'].item(),
            "extinction" : "nan",
            "teff" : main_query['Fe_H_Teff'].item(),
            "logg" : main_query['Fe_H_log_g'].item(),
            "feh" : main_query['Fe_H_Fe_H'].item(),
            "visibility" : months_visible,
            "max_alt" : (max_alt, max_alt_date),
            "exp_time" : exposure,
        }
    )
    return full_data


# ############# ITERATE THROUGH LIST OF STARS
total = len(inputfiledata.values)
colorstr = '\033[1;33m' # for the progress str
ending = "\r" # for the progress str
# generate the dates for the altitude calculations
dates = generate_datetimes()

# iterate through the list
for line in inputfiledata.itertuples():
    name = line[1] # designation
    pulleddata = retrieve_data(name, dates)
    outputfiledata = pd.concat([outputfiledata, pulleddata.to_frame().T], ignore_index=True)
    
    # print the progress bar
    if (line[0] + 1) == total:
        colorstr = '\033[1;32m'
        ending="\n"
    print('Progress: ' + f'{colorstr}' + f'{((line[0]+1)/total)*100:.1f}' +'\033[0m' + ' %',
        sep='',
        end=ending,
        flush=True
    )


# ############# WRITE OUTPUT TO FILE
filename = "star_data.csv"
outputfiledata.to_csv(
    filename,
    header=True,
    index=False,
    sep=",",
)


