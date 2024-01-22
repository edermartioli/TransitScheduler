"""
    Created on Jan 22 2024
    
    Description: This routine schedule exoplanet transits that are observable for a given object
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratório Nacional de Astrofísica, Brasil.

    Simple usage examples:
    
    python transit_scheduler_obj.py --object="HATS-29 b" --start_date="2024-07-01T20:00:00" --end_date="2024-07-08T10:00:00"
    python transit_scheduler_obj.py --object="AU Mic b" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59"
    
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
parser.add_option("-p", "--params", dest="params", help='Parameters yaml file',type='string',default="params.yaml")
parser.add_option("-e", "--exoplanet_catalog", dest="exoplanet_catalog", help='Exoplanet.eu catalog csv file',type='string',default="")
parser.add_option("-i", "--object", dest="object", help='Object ID (e.g. "HATS-29 b")',type='string', default="")
parser.add_option("-1", "--start_date", dest="start_date", help='Start date (ISOT)',type='string', default="2024-04-01T00:00:00")
parser.add_option("-2", "--end_date", dest="end_date", help='End date (ISOT)',type='string', default="2024-09-02T00:00:00")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with transit_scheduler_obj.py -h "); sys.exit(1);

if options.verbose:
    print('Object ID: ', options.object)
    print('Start date: ', options.start_date)
    print('End date: ', options.end_date)
    print('Exoplanet catalog: ', options.exoplanet_catalog)
    print('Output: ', options.output)

# Load parameters from parameters file
params = init_params(options.params)

# Set catalog path
if options.exoplanet_catalog != "":
    params["EXOPLANET_CATALOG"] = options.exoplanet_catalog

# load catalog
exoplanets = ascii.read(params["EXOPLANET_CATALOG"])
#print(exoplanets.columns)

# set date range for transit events
if options.start_date != "" :
    params["START_DATE"] = options.start_date
if options.end_date != "":
    params["END_DATE"] = options.end_date
    
# set object id from command line input
if options.object != "":
    obj_id = options.object
else :
    print("Input a valid exoplanet ID, exiting ... ")
    exit()

# filter exoplanet table to get only the selectected exoplanet
object_info = exoplanets[exoplanets["name"] == obj_id]

# initialize an output table
tbl = Table()

if len(object_info) :
    tbl = add_observable_transits(params, object_info, tbl=tbl, verbose=options.verbose)
else :
    print ("Object ID: {} not found in the database: {}".format(obj_id, options.exoplanet_catalog))

if len(tbl) == 0 :
    print("There are no full observable transits within the selected time range. Exiting ...")
    exit()
        
tbl.sort("TRANSIT_CEN_JD")

print(tbl)

if options.output != "" :
    tbl.write(options.output, overwrite=True)
