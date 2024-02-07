from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.table import Table, vstack
import numpy as np
from scipy import constants

from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import astropy.units as u

import yaml


def add_observable_transits(params, exoplanets_tbl, tbl=Table(), planet_index=0, verbose=False) :

    obj_id = exoplanets_tbl['name'][planet_index]
    per = exoplanets_tbl['orbital_period'][planet_index]
    ra = exoplanets_tbl['ra'][planet_index]
    dec = exoplanets_tbl['dec'][planet_index]
    
    # make sure epoch is defined
    if np.isfinite(exoplanets_tbl['tconj'][planet_index]) :
        ep = exoplanets_tbl['tconj'][planet_index]
    elif np.isfinite(exoplanets_tbl['tzero_tr'][planet_index]) :
        ep = exoplanets_tbl['tzero_tr'][planet_index]
    else :
        print("WARNING: input catalog does not provide a time of transit for object: {}. Exiting ... ".format(obj_id))
        return tbl
    
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
        
    if verbose :
        print("Searching observable transits for object {}/{}: {} RA={} DEC={} P={} d  Gmag={:.2f} Teff={:.0f}".format(planet_index+1, len(exoplanets_tbl), obj_id, ra, dec, per, Gmag, GTeff))

    # add more parameters
    inc = exoplanets_tbl['inclination'][planet_index]
    rp = exoplanets_tbl['radius'][planet_index]
    mstar = exoplanets_tbl['star_mass'][planet_index]
    rstar = exoplanets_tbl['star_radius'][planet_index]
    if not np.isfinite(rstar) :
        rstar = 1.0
    if not np.isfinite(mstar) :
        mstar = 1.0
    ecc = 0
    if np.isfinite(exoplanets_tbl['eccentricity'][planet_index]):
        ecc = exoplanets_tbl['eccentricity'][planet_index]
    om = 90
    if np.isfinite(exoplanets_tbl['omega'][planet_index]):
        om = exoplanets_tbl['omega'][planet_index]
       
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
                                 start_date = params["START_DATE"],
                                 end_date = params["END_DATE"],
                                 date_format='isot',
                                 TimeZone = params["TIMEZONE"],
                                 longitude = params["LONGITUDE"],
                                 latitude = params["LATITUDE"],
                                 altitude = params["ALTITUDE"],
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
        if params["INCLUDE_GAIA_DR3_INFO"] :
            result["GAIA_ID"] = np.full_like(result["FULLY_OBSERVABLE"],GID)
            result["G_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],Gmag)
            result["GAIA_TEFF"] = np.full_like(result["FULLY_OBSERVABLE"],GTeff)
        result["ORB_PERIOD"] = np.full_like(result["FULLY_OBSERVABLE"],per)
        result["TRANSIT_DEPTH"] = np.full_like(result["FULLY_OBSERVABLE"],depth)
        result["TRANSIT_DURATION"] = np.full_like(result["FULLY_OBSERVABLE"],dur)
        result["STAR_MASS"] = np.full_like(result["FULLY_OBSERVABLE"],mstar)
        result["V_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],exoplanets_tbl['mag_v'][planet_index])
        result["TEFF"] = np.full_like(result["FULLY_OBSERVABLE"],exoplanets_tbl['star_teff'][planet_index])

        if verbose :
            print("\t {} full observable transits!".format(len(result)))
            
        # add observable events to the output table
        tbl = vstack([tbl, result])

    return tbl


def add_toi_observable_transits(params, tois_tbl, tbl=Table(), toi_index=0, verbose=False) :
                            
    tic_id = "TIC {}".format(tois_tbl['TIC'][toi_index])
    full_toi_id = "TOI-{}".format(tois_tbl['Full TOI ID'][toi_index])
    toi_id = "TOI-{:.0f}".format(np.floor(tois_tbl['Full TOI ID'][toi_index]))
    
    per = tois_tbl['Orbital Period Value'][toi_index]
    ra = tois_tbl['TIC Right Ascension'][toi_index]
    dec = tois_tbl['TIC Declination'][toi_index]
    
    ep = tois_tbl['Epoch Value'][toi_index] + 2457000.0
    ep_err = tois_tbl['Epoch Error'][toi_index]
    
    tmag = tois_tbl['TMag Value'][toi_index]
    
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
        
    if verbose :
        print("Searching observable transits for object {}/{}: {} {} RA={} DEC={} P={} d  Tmag={:.2f} GTeff={:.0f}".format(toi_index+1, len(tois_tbl), full_toi_id, tic_id, ra, dec, per, tmag, GTeff))

    # calculate transit duration
    dur = tois_tbl['Transit Duration Value'][toi_index]
    dur_err = tois_tbl['Transit Duration Error'][toi_index]
        
    # get transit depth
    depth = tois_tbl['Transit Depth Value'][toi_index] / 1e6
    depth_err = tois_tbl['Transit Depth Error'][toi_index] / 1e6

    # run routine to search for observable transits
    result = observable_transits(ep,
                                 per,
                                 ra,
                                 dec,
                                 dur,
                                 start_date = params["START_DATE"],
                                 end_date = params["END_DATE"],
                                 date_format='isot',
                                 TimeZone = params["TIMEZONE"],
                                 longitude = params["LONGITUDE"],
                                 latitude = params["LATITUDE"],
                                 altitude = params["ALTITUDE"],
                                 minimum_baseline = params["MIN_BASELINE"],
                                 obj_id = tic_id,
                                 verbose = False,
                                 full_transits_only = params["FULL_TRANSITS_ONLY"],
                                 ra_in_hourangle = False )
                          
    # filter by transit depths
    if len(result) and depth > params["MIN_TRANSIT_DEPTH"] and depth < params["MAX_TRANSIT_DEPTH"]:
        
        # Feed information into the results table
        #result["TOI_ID"] = np.full_like(result["FULLY_OBSERVABLE"],toi_id,dtype="str")
        #result["TIC_ID"] = np.full_like(result["FULLY_OBSERVABLE"],tic_id,dtype="str")
        result["TIC_ID"] = np.full_like(result["FULLY_OBSERVABLE"],tois_tbl['TIC'][toi_index])
        result["TOI_ID"] = np.full_like(result["FULLY_OBSERVABLE"],toi_id,dtype='str')
        result["FULL_TOI_ID"] = np.full_like(result["FULLY_OBSERVABLE"],tois_tbl['Full TOI ID'][toi_index])
        result["T_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],tmag)
        if params["INCLUDE_GAIA_DR3_INFO"] :
            result["GAIA_ID"] = np.full_like(result["FULLY_OBSERVABLE"],GID)
            result["G_MAG"] = np.full_like(result["FULLY_OBSERVABLE"],Gmag)
            result["GAIA_TEFF"] = np.full_like(result["FULLY_OBSERVABLE"],GTeff)
        result["ORB_PERIOD"] = np.full_like(result["FULLY_OBSERVABLE"],per)
        result["TRANSIT_DEPTH"] = np.full_like(result["FULLY_OBSERVABLE"],depth)
        result["TRANSIT_DURATION"] = np.full_like(result["FULLY_OBSERVABLE"],dur)

        if verbose :
            print("\t {} full observable transits!".format(len(result)))
            
        # add observable events to the output table
        tbl = vstack([tbl, result])

    return tbl


def init_params(params_filename) :

    # Load parameters from parameters file
    with open(params_filename, 'r') as f:
        params = yaml.safe_load(f)

    return params
    

def reset_params_for_object_search(params) :

    params['MIN_TRANSIT_DEPTH'] = 0.0
    params['MAX_TRANSIT_DEPTH'] = 1.0

    # Define minimum and maximum orbital period (days)
    params['MIN_PERIOD'] = 0.0
    params['MAX_PERIOD'] = 1e10

    # Define minimum and maximum planet mass (MJup)
    params['MIN_MASS'] = 0.0
    params['MAX_MASS'] = 1e10

    # Define minimum and maximum V magnitude
    params['MIN_VMAG'] = -1e10
    params['MAX_VMAG'] = 1e10

    # Define minimum and maximum T magnitude (For TOI objects)
    params['MIN_TMAG'] = -1e10
    params['MAX_TMAG'] = 1e10

    # Define minimum and maximum declination (deg)
    params['MIN_DECLINATION'] = -90
    params['MAX_DECLINATION'] = +90

    # Define minimum and maximum right ascention (deg)
    params['MIN_RA'] = 0
    params['MAX_RA'] = 360
    
    return params
    
def get_sun_alt(observing_time, observer_location) :

    # Calculate the Sun's position at the given observing time
    sun_position = get_sun(observing_time)

    # transform sun position to local coordinates
    loc_sun_position = sun_position.transform_to(AltAz(obstime=observing_time, location=observer_location))
    
    return loc_sun_position.alt
 

def observable_transits(ep, per, ra, dec, dur, start_date="2023-04-01T00:00:00", end_date="2023-12-01T00:00:00", date_format='isot', equinox="J2000", evening_twillight=18.0, morning_twillight=6.0, TimeZone=-3, longitude=-45.5825, latitude=-22.53444, altitude=1864, minimum_baseline=0.0, obj_id="UNKNOWN", full_transits_only=False, ra_in_hourangle=True, verbose=False) :

    # Define reference epoch of conjunction
    ep_time = Time(ep, format='jd', scale='utc')

    if verbose :
        print("-----------------")
        print("INPUT PARAMETERS:")
        print("-----------------")
        print("Object:",obj_id)
        print("Reference epoch of conjunction: JD {} ; ISOT {} ".format(ep_time.jd, ep_time.isot))
        print("Orbital period: {} days".format(per))
        print("Coordinates: RA={} deg DEC={} deg".format(ra,dec))
        print("Transit duration: {} hours".format(dur))
        print("Time range: {} -- {}".format(start_date,end_date))
        print("Observatory location:  East-Longitude={} deg North-Latitude={} deg Elevation={} m".format(longitude,latitude,altitude))
        print("Observatory time zone:  UT + ({} h)".format(TimeZone))
        print("***********************************")

    # defube local time zone and observatory elevation
    locTimeZone = TimeZone*u.hour
    altitude = altitude*u.m
    
    # Define observatory location
    observatory = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude)
    # Define source
    
    ra_unit = u.deg
    if ra_in_hourangle :
        ra_unit = u.hourangle
    
    # set source coordinates
    source = SkyCoord(ra, dec, unit=(ra_unit, u.deg), frame='icrs', equinox=equinox)
    
    # Define initial time to be considered
    t1 = Time(start_date, format=date_format, scale='utc')
    # Define final time to be considered
    t2 = Time(end_date, format=date_format, scale='utc')
    
    # Define maximum number of epochs that fits within the selected time range
    neps = int(np.ceil((t2.jd - t1.jd)/per)) + 2
    
    # calculate time difference between reference transit mid-time and initial date
    delta_t = np.abs(t1.jd - ep_time.jd)

    # Figure out first epoch within selected time range
    first_epoch = ep_time.jd
    
    if ep_time.jd <= t1.jd :
        ep_offset = int(np.floor(delta_t/per))
        first_epoch += (ep_offset - 1) * per
    elif ep_time.jd > t1.jd :
        ep_offset = int(np.ceil(delta_t/per)) + 1
        first_epoch -= ep_offset * per
 
    #print("**** T1: {:.3f} *****".format(t1.jd))
    #print("T0 = {:.2f}".format(first_epoch))
    # create array of transit mid times
    tcs = np.array([])
    for j in range(neps) :
        ttransit = first_epoch + j*per
        #print("Transit {}/{} : {:.3f}".format(j+1, neps, ttransit))
        # make sure the transit falls within the selected time range
        if (ttransit >= t1.jd) and (ttransit <= t2.jd) :
            tcs = np.append(tcs, ttransit)
    #for j in range(len(tcs)) :
    #    print(j, tcs[j])
    #print("**** T2: {:.3f} *****".format(t2.jd))

    # set out-of-transit baseline in units of days
    baseline = TimeDelta(u.day*float(minimum_baseline)/24)

    # set transit duration in units of days
    transit_duration = TimeDelta(u.day*(dur/24))

    # Define the altitude for morning and evening twilight (usually -6 degrees for civil twilight)
    twilight_altitude = -6 * u.deg

    # initialize arrays to store data
    t_ingress, am_ingress = np.array([]), np.array([])
    t_center, t_center_jd, am_center = np.array([]), np.array([]), np.array([])
    t_egress, am_egress = np.array([]), np.array([])
    fully_observable = np.array([])
    obj_ids = np.array([])

    stringfile = ""
    for j in range(len(tcs)) :
    
        tt = Time(tcs[j], format='jd',scale='utc')
        tin = tt - transit_duration/2 - baseline
        teg = tt + transit_duration/2 + baseline
        
        sun_alt_tin = get_sun_alt(tin, observatory)
        sun_alt_teg = get_sun_alt(teg, observatory)
        
        tt_AltAz = AltAz(obstime=tt,location=observatory)
        ti_AltAz = AltAz(obstime=tin,location=observatory)
        te_AltAz = AltAz(obstime=teg,location=observatory)
        
        airmass_tt = source.transform_to(tt_AltAz).secz
        airmass_ti = source.transform_to(ti_AltAz).secz
        airmass_te = source.transform_to(te_AltAz).secz
        
        lt_ingress = tin+locTimeZone
        lt_tt = tt+locTimeZone
        lt_egress = teg+locTimeZone
        
        ingress_h = lt_ingress.datetime.time().hour + lt_ingress.datetime.time().minute / 60 + lt_ingress.datetime.time().second / 3600
        egress_h = lt_egress.datetime.time().hour + lt_egress.datetime.time().minute / 60 + lt_egress.datetime.time().second / 3600
        
        full_transit = False
        if (1.0 < airmass_ti < 2.0 and 1.0 < airmass_te < 2.0) and (sun_alt_tin < twilight_altitude) and (sun_alt_teg < twilight_altitude) :
            full_transit = True
            
        if (1.0 < airmass_tt < 2.0 or 1.0 < airmass_ti < 2.0 or 1.0 < airmass_te < 2.0) and (sun_alt_tin < twilight_altitude or sun_alt_teg < twilight_altitude) :
            stringline = "{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\n".format(full_transit, obj_id, (lt_ingress).isot, airmass_ti, (lt_tt).isot, airmass_tt, (lt_egress).isot, airmass_te)
            
            stringfile += stringline

            t_ingress = np.append(t_ingress,(lt_ingress).isot)
            am_ingress = np.append(am_ingress,airmass_ti)

            t_center = np.append(t_center,(lt_tt).isot)
            t_center_jd = np.append(t_center_jd,tt.jd)
            am_center = np.append(am_center,airmass_tt)

            t_egress = np.append(t_egress,(lt_egress).isot)
            am_egress = np.append(am_egress,airmass_te)

            fully_observable = np.append(fully_observable, full_transit)
            obj_ids = np.append(obj_ids, obj_id)
            
    if verbose :
        print(stringfile)

    tbl = Table()
    
    tbl["OBJ_ID"] = obj_ids
    tbl["RA"] = np.full_like(obj_ids, ra)
    tbl["DEC"] = np.full_like(obj_ids, dec)
    tbl["TRANSIT_CEN_JD"] = t_center_jd
    tbl["TRANSIT_IN"] = t_ingress
    tbl["AIRMASS_IN"] = am_ingress
    tbl["TRANSIT_CEN"] = t_center
    tbl["AIRMASS_CEN"] = am_center
    tbl["TRANSIT_OUT"] = t_egress
    tbl["AIRMASS_OUT"] = am_egress
    tbl["FULLY_OBSERVABLE"] = fully_observable
    
    if full_transits_only :
        tbl = tbl[tbl["FULLY_OBSERVABLE"]==1]
    
    return tbl



def transit_duration(per, inc, ecc, om, rp, mstar, rstar) :
    """
        Calculate the transit duration
        Args:
        per (float): orbital period (day)
        inc (float): orbital inclination (degree)
        ecc (float): eccentricity
        om (float): argument of periastron (degree)
        rp (float): planet radius (Rjup)
        mstar (float): stellar mass (Msun)
        rstar (float): stellar radius (Rsun)
        Returns:
        tdur (float): transit duration (days)
        """
    rsun = 696.34e6
    rs = rstar * rsun
    rjup = 69.911e6
    rp_over_rs = rp * rjup / rs
    
    sma_over_rs = semi_major_axis(per, mstar, rstar, over_rs=True)

    ww = om * np.pi / 180
    ii = inc * np.pi / 180
    ee = ecc
    aa = sma_over_rs
    ro_pt = (1 - ee ** 2) / (1 + ee * np.sin(ww))
    b_pt = aa * ro_pt * np.cos(ii)
    if b_pt > 1:
        b_pt = 0.5
    s_ps = 1.0 + rp_over_rs
    df = np.arcsin(np.sqrt((s_ps ** 2 - b_pt ** 2) / ((aa ** 2) * (ro_pt ** 2) - b_pt ** 2)))
    tdur = (per * (ro_pt ** 2)) / (np.pi * np.sqrt(1 - ee ** 2)) * df
    
    return tdur

def semi_major_axis(per, mstar, rstar, over_rs=False) :
    """
        Calculate the semi-major axis
        Args:
        per (float): orbital period in days
        mstar (float): stellar mass in Msun
        rstar (float): stellar radius in Rsun
        over_rs (bool): switch to output units in relative units to star radius
        Returns:
        a (float): semi-major axis in units of au
        """
        
    G = constants.G
    msun = 1.989e+30
    rsun = 696.34e6
    ms = mstar * msun
    rs = rstar * rsun
    d2s = 24.*60.*60.
    
    a = (((d2s*per * d2s*per * G * ms)/(4. * np.pi * np.pi))**(1/3))
    
    if over_rs :
        return a / rs
    else :
        au_m = 1.496e+11
        return a / au_m
 

def transit_depth(rstar, rp) :
    """
        Calculate the transit depth from planet and star radii
        Args:
        rstar (float): stellar radius in Rsun
        rp (float): planet radius in RJup
        Returns:
        depth (float): transit depth
        """
    rsun = 696.34e6
    rjup = 69.911e6
    
    rp_over_rs = (rp * rjup) / (rstar * rsun)
    depth = rp_over_rs * rp_over_rs
    return depth
    

#########################
#### INPUT TRANSIT PARAMS
#########################
"""
obj_id = '270355392'
tessbjd_zero = 2457000
ep, per = 1411.3951504453178+tessbjd_zero, 5.025651951825772
dur = 2.0
ra, dec = 32.2181,4.3064

start_date="2024-04-01T00:00:00"
end_date="2024-09-02T00:00:00"

result = observable_transits(ep, per, ra, dec, dur, start_date=start_date, end_date=end_date, obj_id=obj_id, verbose=True, full_transits_only=True)



#########################
#### INPUT TRANSIT PARAMS
#########################
obj_id = 'HATS-29'
ep, per = 2457031.95618, 4.6058749
dur = 3.2
ra, dec = "19 00 23.1572835384","-54 53 35.449921500"
#########################

start_date="2024-04-01T00:00:00"
end_date="2024-09-02T00:00:00"

opd_result = observable_transits(ep, per, ra, dec, dur, start_date=start_date, end_date=end_date, obj_id=obj_id, verbose=True, full_transits_only=True)

gemini_south_result = observable_transits(ep, per, ra, dec, dur, start_date=start_date, end_date=end_date, TimeZone=-4, longitude=289.2767, latitude=-30.24166667, altitude=2737, obj_id=obj_id, verbose=True, full_transits_only=True)

print("#########################")
print("OPD: ")
print("#########################")

print(opd_result)


print("#########################")
print("Gemini South: ")
print("#########################")

print(gemini_south_result)
"""
