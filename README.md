# TransitScheduler
Toolkit to schedule transits of exoplanets

Simple usage:
--
## To search all planets in database
```
python transit_scheduler.py --start_date="2024-05-03T18:00:00" --end_date="2024-06-02T07:00:00"  --params="/Users/eder/ExoplanetScheduler/params.yaml"
```
## To search an individual planet
```
python transit_scheduler_obj.py --object="AU Mic b" --start_date="2024-01-01T00:00:00" --end_date="2024-12-31T23:59:59"
```
