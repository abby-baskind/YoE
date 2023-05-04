# -----------------------------------------------------------------------------------------|
# READING BUOY DATA TXT FILES FROM URL AND WRITING TO GOOGLE SHEET                         |
# PMEL CO2 Moorings: https://data.pmel.noaa.gov/pmel/erddap/files/all_pmel_co2_moorings/   |
# Abby Baskind                                                                             |
# 4 May 2023                                                                               |
# -----------------------------------------------------------------------------------------|

from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import xarray as xr
import pandas as pd
import scipy
from datetime import datetime, timedelta
import time
import seaborn
import matplotlib.dates as mdates
import bottleneck as bn
import PyCO2SYS as pyco2
import gsw
import math
import netCDF4 as nc
import requests

from importlib import reload
import warnings
warnings.filterwarnings('ignore')

# -----------------------------------------------------------------------------------------|
# # Set up API credentials
# -----------------------------------------------------------------------------------------|
import gspread
from df2gspread import df2gspread as d2g
from oauth2client.service_account import ServiceAccountCredentials

scope = [
   'https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']
#Name of our Service Account Key
google_key_file = '../phyto/servicecredentials.json'  # Credentials file in phyto folder
credentials = ServiceAccountCredentials.from_json_keyfile_name(google_key_file, scope)
gc = gspread.authorize(credentials)

# -----------------------------------------------------------------------------------------|
# # URL retrieval
# -----------------------------------------------------------------------------------------|
# Individual datasets
datalist = ['ALAWAI_w_meta.txt',
            'BOBOA_w_meta.txt',
            'BTM_w_meta.txt',
            'CAPEARAGO_w_meta.txt',
            'CAPEELIZABETH_w_meta.txt',
            'CCE1_w_meta.txt',
            'CCE2_w_meta.txt',
            'CHABA_w_meta.txt',
            'CHEECAROCKS_w_meta.txt',
            'CHUUK_w_meta.txt',
            'COASTALLA_w_meta.txt',
            'COASTALMS_w_meta.txt',
            'CRESCENTREEF_w_meta.txt',
            'CRIMP1_w_meta.txt',
            'CRIMP2_w_meta.txt',
            'DABOB_w_meta.txt',
            'FIRSTLANDING_w_meta.txt',
            'GAKOA_w_meta.txt',
            'GRAYSREEF_w_meta.txt',
            'GULFOFMAINE_w_meta.txt',
            'HOGREEF_w_meta.txt',
            'ICELAND_w_meta.txt',
            'JKEO_w_meta.txt',
            'KANEOHE_w_meta.txt',
            'KEO_w_meta.txt',
            'KILONALU_w_meta.txt',
            'KODIAK_w_meta.txt',
            'LAPARGUERA_w_meta.txt',
            'M2_w_meta.txt',
            'NH10_w_meta.txt',
            'PAPA_w_meta.txt',
            'SEAK_w_meta.txt',
            'SOFS_w_meta.txt',
            'STRATUS_w_meta.txt',
            'TAO110W_w_meta.txt',
            'TAO125W_w_meta.txt',
            'TAO140W_w_meta.txt',
            'TAO155W_w_meta.txt',
            'TAO165E_w_meta.txt',
            'TAO170W_w_meta.txt',
            'TAO8S165E_w_meta.txt',
            'TWANOH_w_meta.txt',
            'WHOTS_w_meta.txt']

# Base URL
urlbase = 'https://data.pmel.noaa.gov/pmel/erddap/files/all_pmel_co2_moorings/'

# -----------------------------------------------------------------------------------------|
# # Combine base URL with individual datasets to get each txt file
# -----------------------------------------------------------------------------------------|
for i in datalist:
    if i == 'ALAWAI_w_meta.txt':
        URL = urlbase + i
        # first 115 rows are just notes
        df = pd.read_csv(URL, skiprows=115)
    else:
        URL = urlbase + i 
        # first 115 rows are just notes
        df2 = pd.read_csv(URL, skiprows=115)
        df = pd.concat([df, df2])

# -----------------------------------------------------------------------------------------|
# # Drop NaNs and reindex
# -----------------------------------------------------------------------------------------|

# Using only `dropna()` gives ambiguous truth values since it's basically a mask
# `reset_index()` solves this
df0 = df.dropna(how = 'any')
dff = df0.reset_index(drop=True)

# -----------------------------------------------------------------------------------------|
# # Round coordinates to match model coordinates
# -----------------------------------------------------------------------------------------|

# If the coordinate is rounded up, subtract 0.5 for the midpoint of the grid cell.
# If the coordinate is rounded down, add 0.5 for the midpoint of the grid cell.

for i in dff.index:
    lati = dff['latitude'][i]        # ex: lati = 27.8
    roundedlat = round(lati)         # ex: roundedlat = 28
    if roundedlat > lati:            # ex: if 28 > 27.8
        finallat = roundedlat - 0.5  # ex: finallat = 28 - 0.5 = 27.5
    else:
        finallat = roundedlat + 0.5
    longi = dff['longitude'][i]
    roundedlon = round(longi)
    if roundedlon > longi:
        finallon = roundedlon - 0.5
    else:
        finallon = roundedlon + 0.5
    dff['latitude'][i] = finallat
    dff['longitude'][i] = finallon
    
# -----------------------------------------------------------------------------------------|
# # Solve for fugacity and Omega; add to data                                              |
#                                                                                          |
# INPUTS                                                                                   |
# > `pCO2_sw`: partial pressure of CO2 in seawater [uatm] -- type `4`                      |
#
# > `pH_sw`: seawater pH -- type `3`
#
# > `SSS`: sea surface salinity
#
# > `SST`: sea surface temperature [degC]
# -----------------------------------------------------------------------------------------|

results = pyco2.sys(par1=dff['pCO2_sw'],par2=dff['pH_sw'],par1_type=4,par2_type=3,salinity = dff['SSS'], temperature = dff['SST'])
dff['fCO2'] = results['fCO2']
dff['OmegaAr'] = results['saturation_aragonite']

# -----------------------------------------------------------------------------------------|
# # Write to spreadsheet shared to Wang Lab Drive
# -----------------------------------------------------------------------------------------|

#This is the Worksheet ID
spreadsheet_key = '1jfjQkobvAVB4eIjg6sK3lKis1l-_czaeoK7tGZfg8VE'
#This is the sheet name
wks_name = 'Sheet1'
#We upload the data to our Google Sheet. Setting the row_names to False if you did not want the index to be included
d2g.upload(dff, spreadsheet_key, wks_name, credentials=credentials, row_names=False)




