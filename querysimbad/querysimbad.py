"""
Created: 2023-08-03
Last Edit: 2023-09-28
Author: Caleb Brand
Contact: brand.c.research@gmail.com
"""

import os
import time
import numpy as np
import pandas as pd
import datetime as dt
import astropy.units as u
import multiprocessing as mp
from astropy.time import Time
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import astropy.coordinates as coords
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Latitude

####################

class querysimbad():

    def __init__(self):
        # set up extra simbad query fields
        Simbad.add_votable_fields(
            'flux(V)',
            'flux(B)',
            'fe_h',
            'sptype',
            'id',
        )
        Simbad.ROW_LIMIT = 0
        # set up mount john location
        self.MJO = EarthLocation(
            lat=-43.9853*u.deg,
            lon=170.4641*u.deg,
            height=1029*u.m,
        )
        # we specifically choose 28 degrees in the following since that is 
        # approximately the lowest the McLellan telescope can observe in the south
        self.min_observable = Latitude(28, unit='deg')
        # set colors for formatting
        self.yellowcolorstr = '\033[1;33m'
        self.greencolorstr = '\033[1;32m'
        self.bluecolorstr = '\033[1;34m'
        self.cancelformatting = '\033[0m'


    def empty_simbad_cache(self):
        # self explanatory, emtpies the simbad query cache
        for fileobject in os.scandir(Simbad.cache_location):
            os.remove(fileobject.path)
        return


    def get_hd_names(self, data):
        # retrieve the HD catalogue number if one exists from the 
        # id field of the original query
        nameoptions = data.split(",")
        for nameoption in nameoptions:
            if nameoption.strip().startswith("HD"):
                return nameoption.strip().replace(' ','')
        return np.nan


    def generate_datetimes(self, dayobject=None):
        # generate datetimes for the alt calculations
        # note all dates and times are in UTC
        current_year = dt.datetime.now().date().year
        if dayobject:
            start_date = dt.datetime(
                dayobject.year,
                dayobject.month,
                dayobject.day,
                0,
                0
            )
            end_date = dt.datetime(
                dayobject.year,
                dayobject.month,
                dayobject.day + 1,
                0,
                0
            )
        else:
            start_date = dt.datetime(current_year, 1, 1, 13, 0)
            end_date = dt.datetime(current_year, 12, 31, 13, 0)
        datetimes = []
        current_date = start_date
        while current_date <= end_date:
            datetimes.append(current_date)
            if dayobject:
                current_date += dt.timedelta(minutes=15)
            else:
                current_date += dt.timedelta(days=1)
        return datetimes


    def estimate_exposure(self, magnitude):
        # estimate an exposure length for the given magnitude
        # and on based on the following dictionary of previously
        # used exposure times for different magnitudes
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
        if np.isnan(magnitude):
            return np.nan
        rounded = np.round(magnitude)
        if rounded < 0:
            return 100
        elif rounded > 7:
            return 1800
        exposure = exposure_dictionary[rounded]
        return exposure


    def get_twilight_times(self, dayobject):
        # calculate the times that the sun sets and rises (in UTC)
        timeobject = Time(dayobject)
        sun = coords.get_sun(timeobject)
        altitudes = self.altitude_calc_setup([sun.ra], [sun.dec], dayobject)
        sun_motion = pd.DataFrame(data = np.column_stack((self.datetimes, [alt.dms[0] for alt in altitudes])))
        nighttime = sun_motion[sun_motion[1] < -12]
        sun_set = nighttime.iloc[0][0]
        sun_rise = nighttime.iloc[-1][0]
        return sun_set, sun_rise


    def altitude_calc_setup(self, ra, dec, dayobject=None):
        # get the altitudes of a given target across an optionally
        # given day, or get the altitudes across the year
        radec = SkyCoord(
            ra,
            dec,
            frame="icrs",
            unit=(u.hourangle, u.deg),
        )
        if dayobject:
            self.datetimes = self.generate_datetimes(dayobject)
        else:
            self.datetimes = self.generate_datetimes()
        altitudes = []
        for date in self.datetimes:
            altaz = radec.transform_to(AltAz(obstime=date, location=self.MJO))
            altitudes.append(altaz.alt)
        return altitudes


    def calculate_max_alt(self, ra, dec):
        # set up altitude calculations
        altitudes = self.altitude_calc_setup(ra, dec)
        max_alt, alt_index = (
            np.max(altitudes, axis=0),
            np.argmax(altitudes, axis=0)
        )
        max_alt_date = [self.datetimes[indexthing].date().month for indexthing in alt_index]
        # iterate of the each 'star' column in altitudes to find the months visible
        altitudes = pd.DataFrame(data=np.array(altitudes)*u.deg)
        visibility = []
        for column in altitudes:
            dates_visible = altitudes[column][altitudes[column] > self.min_observable].index.values
            when_visible = []
            for index in dates_visible:
                when_visible.append(self.datetimes[index].date().month)
            when_visible = list(set(when_visible))
            if len(when_visible) == 12:
                when_visible = "throughout"
            visibility.append(when_visible)
        return max_alt, max_alt_date, visibility


    def plot_alt(self, designation, day, legend=True):
        # plot the altitude of a target across a given night
        day_object = dt.datetime.strptime(day, '%Y-%m-%d')
        data = self.retrieve_data(designation)
        altitudes = self.altitude_calc_setup(data["ra"], data["dec"], day_object)
        sun_set, sun_rise = self.get_twilight_times(day_object)
        print("plotting...")
        # set up plt.rcParams
        plt.rcParams.update({
            "xtick.minor.visible" : True,
            "ytick.minor.visible" : True,
            "font.size" : 12,
            "font.family" : "serif",
            "legend.loc" : "upper center",
        })

        fig,ax = plt.subplots(figsize=(8,8), layout="tight")
        ax.plot(
            self.datetimes,
            altitudes,
            linewidth=1,
            label=designation,
        )
        ax.axvline(
            sun_set,
            color='gray',
        )
        ax.axvline(
            sun_rise,
            color='gray',
        )
        ax.axvspan(
            sun_set,
            sun_rise,
            alpha=0.3,
            color="gray",
        )
        ax.axhline(
            self.min_observable.dms[0],
            color='red',
            linestyle='--',
            label='McLellan Telescope Horizon Limit'
        )
        ax.grid(alpha=0.3)
        if legend:
            ax.legend()
        ax.set_ylabel("Altitude above horizon")
        ax.set_xlabel(f"UTC, {day}")
        ax.set_ylim(0,90)
        ax.set_xlim(
            sun_set - dt.timedelta(hours=4.5),
            sun_rise + dt.timedelta(hours=4.5),
            )
        ticks = ax.get_xticks()
        new_tick_labels = self.convert_xticks(ticks)
        ax.set_xticks(
            ticks = ticks,
            labels = new_tick_labels,
            rotation=0,
        )
        plt.savefig("visibility.png", dpi=300)
        plt.show()
        return


    def convert_xticks(self, ticks):
        times = np.array(ticks)
        times = times*86400
        times = [dt.datetime.utcfromtimestamp(time).time().isoformat(timespec='minutes') for time in times]
        return times


    def retrieve_data(self, designation):
        # get all the data from the simbad query
        if type(designation) == str:
            main_query = Simbad.query_object(designation)
        else:
            main_query = Simbad.query_objects(designation)
        main_query = main_query.to_pandas()
        hd_names = main_query["ID"].apply(self.get_hd_names)
        vmag = main_query['FLUX_V']
        exposure = vmag.apply(self.estimate_exposure)
        max_alt, max_alt_date, months_visible = self.calculate_max_alt(
            main_query["RA"],
            main_query["DEC"]
        )
        full_data = pd.DataFrame(
            data = {
                "name" : main_query["MAIN_ID"],
                "hd" : hd_names,
                "entered_name" : designation,
                "ra" : main_query["RA"],
                "dec" : main_query["DEC"],
                "type" : main_query["SP_TYPE"],
                "vmag" : vmag,
                "bmv" : main_query['FLUX_B'] - main_query['FLUX_V'],
                "extinction" : np.nan,
                "teff" : main_query['Fe_H_Teff'],
                "logg" : main_query['Fe_H_log_g'],
                "feh" : main_query['Fe_H_Fe_H'],
                "visibility" : months_visible,
                "max_alt" : max_alt,
                "max_alt_month" : max_alt_date,
                "exp_time" : exposure,
            }
        )
        print("batch complete...")
        return full_data


    def query_star_list(self, inputfilepath, batchsize=100, filename="stardata.csv"):
        # create output dataframe
        self.outputfiledata = []
        with open(inputfilepath, newline="\n") as f:
            inputfiledata = f.read().splitlines()
        self.batchsize = batchsize
        self.iterations = int(np.ceil(len(inputfiledata) / self.batchsize))
        # count the number of cores and use up to 6, to prevent spamming simbad with query requests
        cores = min(mp.cpu_count()-1, 6)
        # create list of batches of designations
        names = []
        for iteration in range(self.iterations):
            names.append(inputfiledata[iteration*self.batchsize: iteration*self.batchsize + self.batchsize])
        # print progress
        print(f"starting " + 
            self.bluecolorstr + f"{len(names)}" + self.cancelformatting +
            " batches using " + 
            self.bluecolorstr + f"{cores}" + self.cancelformatting + " cores...")
        starttime = time.perf_counter()
        # run the data retrieval over multiple cpus
        with mp.Pool(cores) as p:
            outputfiledata = p.map(self.retrieve_data,names)
        # concat all the data into a single dataframe
        outputfiledata = pd.concat(outputfiledata)
        # write output to file
        outputfiledata.to_csv(
            filename,
            header=True,
            index=False,
            sep=",",
        )
        # 'update' on progress
        endtime = time.perf_counter()
        print(self.bluecolorstr + "completed" + self.cancelformatting)
        print(f"Query time: {endtime - starttime:.2f} secs")



    def plot_alt_year(self, designation, legend=True):
        # plot the altitude of a target across the year
        data = self.retrieve_data(designation)
        altitudes = self.altitude_calc_setup(data["ra"], data["dec"])
        print("plotting...")
        # set up plt.rcParams
        plt.rcParams.update({
            "xtick.minor.visible" : True,
            "ytick.minor.visible" : True,
            "font.size" : 16,
            "font.family" : "serif",
            "legend.loc" : "upper center",
        })

        fig,ax = plt.subplots(figsize=(16,8), layout="tight")
        ax.plot(
            self.datetimes,
            altitudes,
            linewidth=1,
            label=designation,
        )
        ax.axhline(
            self.min_observable.dms[0],
            color='red',
            linestyle='--',
            label='McLellan Telescope Horizon Limit'
        )
        ax.grid(alpha=0.3)
        if legend:
            ax.legend()
        ax.set_ylabel("Altitude above horizon")
        ax.set_xlabel("Date")
        ax.set_title("Altitude as calulated at 13:00 UT for each night of the year.")
        ax.set_ylim(0,90)
        ax.set_xlim(self.datetimes[0],self.datetimes[-1])
        plt.savefig("yearroundvisibility.png", dpi=300)
        plt.show()
        return