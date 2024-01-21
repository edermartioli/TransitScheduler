"""
    Created on Jan 5 2024
    
    Description: This routine schedule exoplanet transits that are observable
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica, Brasil.

    Simple usage examples:
    
    python transit_scheduler.py --start_date="2024-05-03T18:00:00" --end_date="2024-06-02T07:00:00"  --params="/Users/eder/ExoplanetScheduler/params.yaml"
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import sys, os

from optparse import OptionParser

from astropy.table import Table, vstack
import numpy as np
from astropy.io import ascii

from obs_transits_lib import observable_transits, transit_duration, transit_depth

from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import astropy.units as u

import yaml

parser = OptionParser()
parser.add_option("-o", "--output", dest="output", help='Output',type='string', default="")
parser.add_option("-p", "--params", dest="params", help='Parameters yaml file',type='string',default="params.yaml")
parser.add_option("-e", "--exoplanet_catalog", dest="exoplanet_catalog", help='Exoplanet.eu catalog csv file',type='string',default="")
parser.add_option("-1", "--start_date", dest="start_date", help='Start date (ISOT)',type='string', default="2024-04-01T00:00:00")
parser.add_option("-2", "--end_date", dest="end_date", help='End date (ISOT)',type='string', default="2024-09-02T00:00:00")

parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with transit_scheduler.py -h "); sys.exit(1);

if options.verbose:
    print('Exoplanet catalog: ', options.exoplanet_catalog)
    print('Output: ', options.output)
    print('Start date: ', options.start_date)
    print('End date: ', options.end_date)

# Load parameters from parameters file
with open(options.params, 'r') as f:
    params = yaml.safe_load(f)

# Set initial site as OPD
TimeZone = params["TIMEZONE"]
longitude = params["LONGITUDE"]
latitude = params["LATITUDE"]
altitude = params["ALTITUDE"]

# Set catalog path
exoplanet_catalog = params["EXOPLANET_CATALOG"]
if options.exoplanet_catalog != "":
    exoplanet_catalog = options.exoplanet_catalog

# load catalog
exoplanets = ascii.read(exoplanet_catalog)
#print(exoplanets.columns)

# set date range for transit events
start_date = params["START_DATE"]
end_date = params["END_DATE"]
if options.start_date != "" :
    start_date = options.start_date
if options.end_date != "":
    end_date = options.end_date


# filter exoplanet database to match search criteria
transiting_planets = exoplanets[(~exoplanets["radius"].mask) &
                                (~exoplanets["mass"].mask) &
                                (~exoplanets["orbital_period"].mask) &
                                (~(exoplanets["detection_type"]=="Imaging")) &
                                (exoplanets["mass"] > params["MIN_MASS"]) &
                                (exoplanets["mass"] < params["MAX_MASS"]) &
                                (exoplanets["orbital_period"] > params["MIN_PERIOD"]) &
                                (exoplanets["orbital_period"] < params["MAX_PERIOD"]) &
                                (exoplanets["dec"] > params["MIN_DECLINATION"]) &
                                (exoplanets["dec"] < params["MAX_DECLINATION"]) &
                                (exoplanets["ra"] > params["MIN_RA"]) &
                                (exoplanets["ra"] < params["MAX_RA"]) &
                                (exoplanets["mag_v"] > params["MIN_VMAG"]) &
                                (exoplanets["mag_v"] < params["MAX_VMAG"])]

# initialize output table
tbl = Table()

# start loop over each selected object in database
for i in range(len(transiting_planets)) :

    obj_id = transiting_planets['name'][i]
    per = transiting_planets['orbital_period'][i]
    ra = transiting_planets['ra'][i]
    dec = transiting_planets['dec'][i]
    
    # make sure epoch is defined
    if np.isfinite(transiting_planets['tconj'][i]) :
        ep = transiting_planets['tconj'][i]
    elif np.isfinite(transiting_planets['tzero_tr'][i]) :
        ep = transiting_planets['tzero_tr'][i]
    else :
        print("Catalog does not provide a time of transit for object: {}. Skipping ... ".format(obj_id))
        continue
    
    GID = -1
    Gmag, GTeff = np.nan, np.nan
    # Below is to retrieve GAIA DR3 information for this object, it may be slow
    if params["INCLUDE_GAIA_DR3_INFO"] :
        try :
            coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
            width, height = u.Quantity(0.01, u.deg), u.Quantity(0.01, u.deg)
            Gaia.ROW_LIMIT = 5
            r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
            GID = r['source_id'][0]
            Gmag, GTeff = r['phot_g_mean_mag'][0], r['teff_val'][0]
        except :
            print("WARNING: Could not retrieve information Gaia DR3 catalog.")
        
    print("Searching observable transits for object {}/{}: {} RA={} DEC={} P={} d  Gmag={:.2f} Teff={:.0f}".format(i+1, len(transiting_planets), obj_id, ra, dec, per, Gmag, GTeff))

    # add more parameters
    inc = transiting_planets['inclination'][i]
    rp = transiting_planets['radius'][i]
    mstar = transiting_planets['star_mass'][i]
    rstar = transiting_planets['star_radius'][i]
    ecc = 0
    if np.isfinite(transiting_planets['eccentricity'][i]):
        ecc = transiting_planets['eccentricity'][i]
    om = 90
    if np.isfinite(transiting_planets['omega'][i]):
        om = transiting_planets['omega'][i]
       
    # calculate transit duration
    dur = transit_duration(per, inc, ecc, om, rp, mstar, rstar) * 24
    if np.isfinite(rp) and not np.isfinite(dur) :
        dur = transit_duration(per, 90, 0, 90, rp, 1, 1) * 24
    elif not np.isfinite(rp) and not np.isfinite(dur):
        dur = transit_duration(per, 90, 0, 90, 0.03, 1, 1) * 24
        
    # run routine to search for observable transits
    result = observable_transits(ep,
                                 per,
                                 ra,
                                 dec,
                                 dur,
                                 start_date = start_date,
                                 end_date = end_date,
                                 date_format='isot',
                                 TimeZone = TimeZone,
                                 longitude = longitude,
                                 latitude = latitude,
                                 altitude = altitude,
                                 minimum_baseline = params["MIN_BASELINE"],
                                 obj_id = obj_id,
                                 verbose = False,
                                 full_transits_only = params["FULL_TRANSITS_ONLY"],
                                 ra_in_hourangle = False )
                          
    # calculate transit depth
    depth = transit_depth(rstar, rp)

    # filter by transit depths
    if len(result) and depth > params["MIN_TRANSIT_DEPTH"] and depth < params["MAX_TRANSIT_DEPTH"]:
        
        # Feed information into the results table
        result["GAIA_ID"] = np.full_like(result["FULLY_OBSERVABLE"],GID)
        result["G_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],Gmag)
        result["GAIA_TEFF"] = np.full_like(result["FULLY_OBSERVABLE"],GTeff)
        result["ORB_PERIOD"] = np.full_like(result["FULLY_OBSERVABLE"],per)
        result["TRANSIT_DEPTH"] = np.full_like(result["FULLY_OBSERVABLE"],depth)
        result["TRANSIT_DURATION"] = np.full_like(result["FULLY_OBSERVABLE"],dur)
        result["STAR_MASS"] = np.full_like(result["FULLY_OBSERVABLE"],mstar)
        result["V_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],transiting_planets['mag_v'][i])
        result["TEFF"] = np.full_like(result["FULLY_OBSERVABLE"],transiting_planets['star_teff'][i])

        # add observable events to the output table
        print("\t {} full observable transits!".format(len(result)))
        tbl = vstack([tbl, result])
        
if len(tbl) == 0 :
    print("There are no full observable transits within the selected time range. Exiting ...")
    exit()
        
# sort final table by JD to make sure the events are ordered by date and not by object
tbl.sort("TRANSIT_CEN_JD")

# print table of observable transit events
if options.verbose :
    print(tbl)

# save table of observable transit events to file
if options.output != "" :
    tbl.write(options.output, overwrite=True)
