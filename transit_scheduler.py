"""
    Created on Jan 22 2024
    
    Description: This routine schedule exoplanet transits that are observable
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica, Brasil.

    Simple usage examples:
    
    python transit_scheduler.py --start_date="2024-05-03T18:00:00" --end_date="2024-06-02T07:00:00"  --params="/Users/eder/ExoplanetScheduler/params.yaml"
    
    python transit_scheduler.py --object="AU Mic b" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59"
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import sys, os

from optparse import OptionParser

from astropy.table import Table
import numpy as np
from astropy.io import ascii

from obs_transits_lib import add_observable_transits, init_params

parser = OptionParser()
parser.add_option("-o", "--output", dest="output", help='Output',type='string', default="")
parser.add_option("-i", "--object", dest="object", help='Object ID (e.g. "HATS-29 b")',type='string', default="")
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
params = init_params(options.params)

# Set catalog path
if options.exoplanet_catalog != "":
    params["EXOPLANET_CATALOG"] = options.exoplanet_catalog

# load exoplanet.eu catalog
exoplanets = ascii.read(params["EXOPLANET_CATALOG"])

# set date range for transit events
if options.start_date != "" :
    params["START_DATE"] = options.start_date
if options.end_date != "":
    params["END_DATE"] = options.end_date

if options.object == "" :
    ### SURVEY MODE ####
    # filter exoplanet database to match search criteria
    selected_planets = exoplanets[(~exoplanets["radius"].mask) &
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

else :
    ### OBJECT MODE ####
    # filter exoplanet table to get only the selectected exoplanet
    selected_planets = exoplanets[exoplanets["name"] == options.object]
    if len(transiting_planets) == 0 :
        print ("Object ID: '{}' not found in the database: {}".format(options.object, params["EXOPLANET_CATALOG"]))
        exit()

# initialize output table
tbl = Table()

# start loop over each selected object in database
for i in range(len(selected_planets)) :
    tbl = add_observable_transits(params, selected_planets, tbl=tbl, planet_index=i, verbose=options.verbose)
    
if len(tbl) == 0 :
    print("There are no observable transits within the selected time range. Exiting ...")
    exit()
        
# sort final table by JD to make sure the events are ordered by date and not by object
tbl.sort("TRANSIT_CEN_JD")

# print table of observable transit events
print(tbl)

# save table of observable transit events to file
if options.output != "" :
    tbl.write(options.output, overwrite=True)
