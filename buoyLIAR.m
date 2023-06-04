clear all

% *************************************************************************
% This script calculates TA for buoy data that has pCO2 but not pH. This
% data will allow us to have sufficient inputs to solve the carbonate
% system for fugacity, pH, and omega. 
% 
% To calculate TA, we use the LIAR method. Data was originally scraped
% using Python and then downloaded as a csv, in order to load into MATLAB
% for LIAR. After calculating TA, we write everything back into a csv with
% TA, which can then be loaded back into the usual Python scripts for
% analysis.
% 
% Abby Baskind, 6/1/2023
% *************************************************************************

% Add path to LIAR.m
addpath '/Users/akbaskind/Documents/MATLAB/LIRs-master';
% Add path to seawater
addpath '/Users/akbaskind/Documents/MATLAB/LIRs-master/seawater_ver3_3.1';
% Add path to csvwrite_with_headers
addpath '/Users/akbaskind/Documents/MATLAB/csvwrite_with_headers';

% Read CSV as a table and extract data from columns
% readcsv() only takes numeric values
df_LIAR = readtable('/Users/akbaskind/Desktop/YoE/buoy_LIAR.csv');
index = df_LIAR.Var1;
datetime = df_LIAR.datetime_utc;
SSS = df_LIAR.SSS;
SST = df_LIAR.SST;
pCO2 = df_LIAR.pCO2_sw;
site = df_LIAR.site;
lat = df_LIAR.latitude;
lon = df_LIAR.longitude;
depth = df_LIAR.depth;

% combine lon, lat, and depth (in that order) for first argument in LIAR
coords = [lon, lat, depth];
% combine SSS and SST for second argument
% Types are 1 and 7
measurements = [SSS, SST];

% Calculate alkalinity using LIAR method
TA = LIAR(coords, measurements, [1,7],'VerboseTF',false);

% df_out is our output data table containing
% Index
% Site Name
% Date Time
% Latitude (not rounded)
% Longitude (not rounded)
% Sea surface salinity
% sea surface temperature
% pCO2
% TA calculated from SSS and SST in mmol/kg
df_out = table(index, site, datetime, lat, lon, SSS, SST, pCO2, TA);

% Rename the columns
df_out = renamevars(df_out,["lat","lon", "pCO2", "TA"],["Latitude","Longitude", "pCO2_sw", "Total Alkalinity [mmol/kg]"]);
% Save to csv
writetable(df_out,'/Users/akbaskind/Desktop/YoE/buoy_LIAR_out.csv','Delimiter',',');
