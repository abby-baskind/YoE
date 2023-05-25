# Model Retrieval 

# `GFDL-ESM4.gr.historical.Omon.r1i1p1f1`

The selected CMIP6 model is GFDL-ESM4 (Dunne et al., 2020; Stock et al., 2020). It has been regridded onto a 1°x1° lat-lon grid. The `historical` experiment on the `Omon` table runs from 1-15-1850 to 12-15-2014 with data provided monthly. I have selected the first ensemble member, `r1i1p1f1`, as that is the only one provided.

Model output was retrieved and processed on Pangeo [(Abernathy et al., 2021, Computing in Science and Engineering)](https://par.nsf.gov/servlets/purl/10287683) with a variety of key open-source Python packages. All these packages are preinstalled on the Pangeo platform.
* `xmip` [(Busecke, 2021)](https://cmip6-preprocessing.readthedocs.io/en/latest/)
* `xesmf` [(Zhuang, 2020)](https://xesmf.readthedocs.io/en/v0.6.3/)
* `xarray` (Hamman & Hoyer, 2017, J. Open Research Software)
* 'PyCO2SYS`(Humphreys et al., 2022, Geosci. Model Dev.)
* `gsw` (McDougall & Barker, 2011, Scor/Iapso WG)


## Variables

I selected the following variables. For 3D variables (i.e. not surface variables), I selected the first depth level, which corresponds to 2.5m depth. In a future version, I should look for surface variables to replace these.
* `talk`: Total Alkalinity [mol m-3]
* `co3sataragos`: Surface Mole Concentration of Carbonate Ion in Equilibrium with Pure Aragonite in Sea Water [mol m-3]
* `fgco2`: Surface Downward Flux of Total CO2 [kg m-2 s-1]
* `co3os`: Surface Carbonate Ion Concentration [mol m-3]
* `sos`: Sea Surface Salinity [0.001]
* `ph`: ph [Total Scale]
* `tos`: Sea Surface Temperature [degC]
* `dissic`: Dissolved Inorganic Carbon Concentration [mol m-3]

### Calculated Variables

Some variables are either in inconvenient units or not provided in the model. To convert from mol m-3 to umol kg-1, I first calculated density as potential density referenced to the surface using the `gsw` package and then used density as part of a conversion factor `conv`:
```python
ds['sigma0'] = gsw.sigma0(ds['sos'],ds['tos'])

conv = 1e6/(1000 + ds['sigma0'])
```
Fugacity and Omega_arag were not provided by the model. I calculated these using DIC, TA, temperature, and salinity, using `PyCO2SYS`. I based this method on [Terhaar et al. (2021, Biogeosciences)](https://bg.copernicus.org/articles/18/2221/2021/bg-18-2221-2021-discussion.html).

```python
results = pyco2.sys(par1=ds['talk']*conv,par2=ds['dissic']*conv,
                       par1_type=1,par2_type=2, salinity = ds['sos'], temperature = ds['tos'])
``` 
## Detrending Data

As of 6 May 2023 I am detrending from 2000

The manuscript suggests detrending "1.89 µatm yr-1 for fCO2, -0.0018 yr-1 for pH and -0.0078 yr-1 for Ω from 1980s to 2010s (Bates et al., 2014)." For the model and for updated SOCAT data, I referenced to 1980, as suggested by the manuscript. The detrending in the model data is done as follows:

```python
ds['fugacity_detrended'] = (ds['fugacity'] -  1.89 * (ds['fugacity'].time.dt.year - 2000)) * xr.ones_like(ds['talk'])
ds['ph_detrended'] = (ds['ph'] + 0.0018 * (ds['ph'].time.dt.year - 2000)) * xr.ones_like(ds['talk'])
ds['omega_detrended'] = (ds['omega'] + 0.0078 * (ds['omega'].time.dt.year - 2000)) * xr.ones_like(ds['talk'])
```
SOCAT detrending done as follows:

```Matlab
REF_YEAR = 2000;
tem2_lon_lat(:,SOCAT.fCO2rec_detrended)=tem2_lon_lat(:,SOCAT.fCO2rec)-1.89*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
tem2_lon_lat(:,SOCAT.pH_detrended)=tem2_lon_lat(:,SOCAT.pH)+0.0018*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
tem2_lon_lat(:,SOCAT.OmegaAr_detrended)=tem2_lon_lat(:,SOCAT.OmegaAr)+0.0078*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
```

## Time Selection
Only after calculating new variables and detrending the data did I take a mean and standard deviation over time. 

### Present: 1-15-1995 to 12-15-2014

I selected the last 20 years of the model. Using the historical simulation.

### Past: 1-15-1850 to 12-15-1879

For a reconstruction of the past, I recommend the first available 20 years of the model. Remember, the historical model starts in 1850. While a longer climatology (e.g. 50 years) might be more robust, given how late the model starts, I contend climatologies past 1880 would begin showing signs of anthropogenic influence. Moreover, time periods longer than 20 years overwhelm the large Pangeo kernel. If necessary, I can upgrade to a huge kernel, but that is a lot of computing power.

### Future

I have yet to decide which model run to use for our future projection. SSP experiments on `Omon` tables start at 1-15-2015 are initialized from the historical experiment, so they form a continuous record, but their forcings are not based on observation.

There are multiple SSPs to choose from. SSPs are first categorized by scenario numbers 1 through 5, with SSP1 being the most sustainable model, SSP2 being a middle of the road model, SSP3 being a model characterized by regional rivalry, SSP4 being a model characterized by widespread inequity, and SSP5 being a model of continued fossil-fuel development. 

The SSP is then further characterized by a 2100 climate forcing similar to the previous RCPs: 2.6 W/m2, 4.5 W/m2, 7.0 W/m2, and 8.5 W/m2. Models are required to publish "Tier 1" experiments, which are SSP5-8.5, SSP3-7.0, SSP2-4.5, and SSP1-2.6. Further description of the SSP experiments can be found in [O'Neill et al. (2016, Geosci. Model Dev.)](https://gmd.copernicus.org/articles/9/3461/2016/).

# SOCAT

## Updated SOCAT retrieval

I updated the SOCAT processing script (MATLAB) written by Hongjie (`VAR_SOCAT_test.m`). The changes largely reflect updates in the CO2SYS MATLAB package. This means we added arguments for NH4 and H2S, both of which are set to 0 in this context. As well, we added arguments for the boron-to-salinity ratio, set to `2` to reference Lee 2010, and for the KHF dissociation constant, set to `2` to reference Perez & Fraga, 1987. The updated CO2SYS uses filler values of -999, so we added in a line that checks for negative values of key CO2SYS outputs and replaced those with `NaN`.

The use of the deprecated CSIRO package in LIAR remains.

## Restructured SOCAT data for Python

I restructured the SOCAT data calculated by Hongjie using `xarray` to make the structure more similar to the model output. As part of this restructuring, I had to adjust latitude and longitude. The model's lat-lon points are half-degrees (i.e. -89.5, -88.5, -87.5, etc.), whereas SOCAT's lat-lon points are whole degrees (i.e. -90, -89, -88, etc.). So I could compare like grid cells, I added 0.5 to SOCAT's coordinates. This may cause a slight error, but I content the error is minimal and otherwise hard to avoid, as the model does not provide data on the exact coordinates as SOCAT.

# [Buoy Data](https://data.pmel.noaa.gov/pmel/erddap/files/all_pmel_co2_moorings/)

Buoy data retrieval is automated by `buoydata_retrieval.py`. Output is written to a Google sheet in the Wang Lab shared drive. fCO2 and OmegaAr were added using PyCO2SYS.



