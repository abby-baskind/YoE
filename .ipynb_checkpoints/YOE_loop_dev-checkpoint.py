from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import xarray as xr
import pandas as pd
import scipy
from datetime import datetime, timedelta
import time
import matplotlib.dates as mdates
import gsw
import math
import netCDF4 as nc
from importlib import reload
import warnings
warnings.filterwarnings('ignore')
import scipy.optimize as opt 

# Open historical dataset
dd = xr.open_mfdataset('/Users/akbaskind/Desktop/YoE/GFDLESM4historical_1850_2000_means.nc')
# Subset to, like, 9 grid cells
dd = dd.sel(x = slice(127,130), y = slice(0,3))
# ---------------------------------------------------------------------------------------------
# OUT: dd
# <xarray.Dataset>
# Dimensions:         (t: 15, y: 3, x: 3)
# Coordinates:
#   * t               (t) int64 1850 1860 1870 1880 1890 ... 1960 1970 1980 1990
#   * y               (y) float64 0.5 1.5 2.5
#   * x               (x) float64 127.5 128.5 129.5
#     lat_bounds      (y, x) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     lon_bounds      (x, y) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     lon             (x, y) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     lat             (x, y) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     lon_verticies   (x, y) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     lat_verticies   (x, y) float64 dask.array<chunksize=(3, 3), meta=np.ndarray>
#     member_id       object ...
#     dcpp_init_year  float64 ...
# Data variables:
#     fCO2            (t, y, x) float64 dask.array<chunksize=(15, 3, 3), meta=np.ndarray>
#     OmegaAr         (t, y, x) float64 dask.array<chunksize=(15, 3, 3), meta=np.ndarray>
#     ph              (t, y, x) float64 dask.array<chunksize=(15, 3, 3), meta=np.ndarray>
# ---------------------------------------------------------------------------------------------

# Create super coarse time array
t = np.arange(1850,2050,25)
# OUT: array([1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025])

# Empty dataset like dd for output data
ddict = xr.Dataset()
ddict = ddict.assign_coords({'y': dd.y, 'x': dd.x, 'time':t})
dat = np.zeros((len(ddict.y), len(ddict.x), len(t)))
dat[:] = np.nan
ddict['dummy'] = xr.Variable(dims = ['y','x','time'], data = dat)
ddict

def quintic_polynomial(px, t): 
    return px[0] + px[1] * t + px[2] * t**2 + px[3] * t**3 + px[4] * t**4 + px[5] * t**5 
    solution = opt.root(quintic_polynomial(px,t), [1,1,1,1,1,1]) 
    return solution.t

# Loop through latitudes i
for i in dd.x:
    # Loop through longitudes j
    for j in dd.y:
        # With the xarray dataset dd select data at lat i and lon j
        ds = dd.sel(x = i, y = j)
        # If data at (i,j) is Nan, dummy var is Nan
        if math.isnan(ds.ph.values.tolist()[0]):
            ddict["dummy"] = np.nan
        # If not nan
        else:
            # Create 5th order polynomial model for the pH data
            px = np.poly1d(np.polyfit(ds.t, ds.ph, 5))
            # Solve fifth order polynomial for pH at times t
            sols = quintic_polynomial(px, t)
            # n: index for pH in dummy variable, corresponding in t
            n = 0
            # Loop through times t
            for k in t:
                # In new dataset, where lat = i, lon = j, and time = k, save the pH for time at index n as "dummy"
                ddict["dummy"] = xr.where((ddict.coords["x"] == i) & (ddict.coords["y"] == j) & (ddict.coords["time"] == k), sols[n], ddict["dummy"])
                # Increment n
                n += 1