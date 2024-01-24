# TransitScheduler
Toolkit to schedule transits of exoplanets.

Required libraries:
```
numpy version='1.24.3'
scipy version='1.11.3'
astropy version='5.3.4'
astroquery version='0.4.6'
yaml version='6.0.1'
optparse version='1.5.3'
```

Simple usage:
--
To search all planets in exoplanet.eu database:
```
python transit_scheduler.py --start_date="2024-05-03T18:00:00" --end_date="2024-06-02T07:00:00"
```
To search an specific planet in exoplanet.eu database:
```
python transit_scheduler.py --object="AU Mic b" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59"
```
To search all exoplanet candidates in TOI database:
```
python transit_scheduler_toi.py --start_date="2024-05-03T18:00:00" --end_date="2024-05-06T07:00:00" --output="observable_TOIs_in_may.csv" -v
```
To search an specific planet in TOI database:
```
python transit_scheduler_toi.py --toi="TOI-1066.01" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59" -v
```

