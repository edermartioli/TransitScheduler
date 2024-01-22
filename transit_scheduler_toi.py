"""
    Created on Jan 22 2024
    
    Description: This routine schedule TOI transits that are observable
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica, Brasil.

    Simple usage examples:
    
    python transit_scheduler_toi.py --start_date="2024-05-03T18:00:00" --end_date="2024-06-02T07:00:00"  --params="/Users/eder/ExoplanetScheduler/params.yaml"
    
    python transit_scheduler_toi.py --toi="TOI-3568.01" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59"
    
    python transit_scheduler_toi.py --start_date="2024-05-03T18:00:00" --end_date="2024-05-06T07:00:00" --output="observable_TOIs_in_may.csv" -v
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

from obs_transits_lib import add_toi_observable_transits, init_params, reset_params_for_object_search

parser = OptionParser()
parser.add_option("-o", "--output", dest="output", help='Output',type='string', default="")
parser.add_option("-i", "--toi", dest="toi", help='Object ID (e.g. "TOI-3568.01")',type='string', default="")
parser.add_option("-p", "--params", dest="params", help='Parameters yaml file',type='string',default="params.yaml")
parser.add_option("-t", "--toi_catalog", dest="toi_catalog", help='TOI catalog csv file',type='string',default="")
parser.add_option("-1", "--start_date", dest="start_date", help='Start date (ISOT)',type='string', default="2024-04-01T00:00:00")
parser.add_option("-2", "--end_date", dest="end_date", help='End date (ISOT)',type='string', default="2024-09-02T00:00:00")

parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with transit_scheduler.py -h "); sys.exit(1);

if options.verbose and options.output != "":
    print('Output: ', options.output)

# Load parameters from parameters file
params = init_params(options.params)

# Set catalog path
if options.toi_catalog != "":
    params["TOI_CATALOG"] = options.toi_catalog

if options.verbose:
    print('TOI catalog: ', params["TOI_CATALOG"])
    
# load exoplanet.eu catalog
tois = ascii.read(params["TOI_CATALOG"])

# set date range for transit events
if options.start_date != "" :
    params["START_DATE"] = options.start_date
if options.end_date != "":
    params["END_DATE"] = options.end_date
    
if options.verbose:
    print('Start date: ', params["START_DATE"])
    print('End date: ', params["END_DATE"])

# initialize output table
tbl = Table()

if options.toi == "" :
    ### SURVEY MODE ####
    # filter exoplanet database to match search criteria
    selected_planets = tois[(tois["Orbital Period Value"] > params["MIN_PERIOD"]) &
                            (tois["Orbital Period Value"] < params["MAX_PERIOD"]) &
                            (tois["TIC Declination"] > params["MIN_DECLINATION"]) &
                            (tois["TIC Declination"] < params["MAX_DECLINATION"]) &
                            (tois["TIC Right Ascension"] > params["MIN_RA"]) &
                            (tois["TIC Right Ascension"] < params["MAX_RA"]) &
                            (tois["TMag Value"] > params["MIN_TMAG"]) &
                            (tois["TMag Value"] < params["MAX_TMAG"])]
else :
    params = reset_params_for_object_search(params)
    toi_number = float(options.toi.replace("TOI-",""))
    ### OBJECT MODE ####
    # filter exoplanet table to get only the selectected exoplanet
    selected_planets = tois[tois['Full TOI ID'] == toi_number]
    if len(selected_planets) == 0 :
        print ("'{}' not found in the database: {}".format(options.toi, params["TOI_CATALOG"]))
        exit()

# start loop over each selected object in database
for i in range(len(selected_planets)) :
    tbl = add_toi_observable_transits(params, selected_planets, tbl=tbl, toi_index=i, verbose=options.verbose)

if len(tbl) == 0 :
    print("There are no TOI observable transits within the selected time range. Exiting ...")
    exit()
        
# sort final table by JD to make sure the events are ordered by date and not by object
tbl.sort("TRANSIT_CEN_JD")

# print table of observable transit events
print(tbl)

# save table of observable transit events to file
if options.output != "" :
    tbl.write(options.output, overwrite=True)
