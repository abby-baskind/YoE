clear all;

% Add path to LIAR.m
addpath '/Users/akbaskind/Documents/MATLAB/LIRs-master';
% Add path to seawater
addpath '/Users/akbaskind/Documents/MATLAB/LIRs-master/seawater_ver3_3.1';
% Add path to CO2SYS
addpath '/Users/akbaskind/Documents/MATLAB/github_repo';
addpath '/Users/akbaskind/Documents/MATLAB/github_repo/comparisons';
addpath '/Users/akbaskind/Documents/MATLAB/github_repo/comparisons/data';
addpath '/Users/akbaskind/Documents/MATLAB/github_repo/comparisons/derivatives';
addpath '/Users/akbaskind/Documents/MATLAB/github_repo/main';
addpath '/Users/akbaskind/Documents/MATLAB/github_repo/v2_0_5_compatible';
% Add path to csvwrite_with_headers
addpath '/Users/akbaskind/Documents/MATLAB/csvwrite_with_headers';
% Add path to SOCAT
addpath '/Users/akbaskind/Documents/MATLAB/SOCAT';


par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par2type =    5; % The first parameter supplied is of type "5", which is "fCO2"
sil      =    0; % Concentration of silicate  in the sample (in umol/kg)
po4      =    0; % Concentration of phosphate in the sample (in umol/kg)
nh4      =    0; % Concentration of NH4 in the sample (in umol/kg)
h2s      =    0; % Concentration of H2S in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    13;% Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % 1 = KSO4 of Dickson   (PREFERRED)
kfc      =    2; % 2 = KF of Perez & Fraga, 1987  (PREFERRED)
tb       =    2; % TB of Lee 2010  (PREFERRED)

% Load selected data fields from SOCAT
load SOCATV2022.mat;
dat=[longitude,latitude,yr,mon,day,hh,mm,ss,SST,sal,fCO2rec];

% Create intermediate data repository with following headers
headers={'lon','lat', 'yr' 'mon','day','hh','mm','ss','SST','SSS','fCO2rec','TA_Liar','OmegaAr','OmegaCa','pH','Hfree','fCO2rec_detrended','OmegaAr_detrended','pH_detrended', 'dist_to_land','fCO2rec_detrended_min','OmegaAr_detrended_min','pH_detrended_min','fCO2rec_detrended_max','OmegaAr_detrended_max','pH_detrended_max'};
 for n=1:size(headers,2);
    SOCAT.(headers{1,n})=n;
 end  
 
 % repository for final data
 % final data are results from SOCAT
 VAR_SOCAT=[];
 
 k=0

% loop thorugh longitudes
for lon=0:1:360;
% for lon=0:1:5;    
        % select 1 longitude grid
        loc1=find((dat(:,SOCAT.lon))<lon+1 & (dat(:,SOCAT.lon))>=lon);
        tem1_lon=dat(loc1,:);               % tem1_lon: all data with given lon
        
        % loop through latitudes
        for lat=-90:1:90;
%         for lat=40:1:42;
            % select 1 latitude grid cell from pre-selected longitudes
            loc2=find((tem1_lon(:,SOCAT.lat))<lat+1 & (tem1_lon(:,SOCAT.lat))>=lat);
            tem2_lon_lat=tem1_lon(loc2,:);  % tem2_lon_lat: all data with given lon and lat
%             tem2_lon_lat(:,12) = datenum([tem2_lon_lat(:,SOCAT.yr),tem2_lon_lat(:,SOCAT.mon),tem2_lon_lat(:,SOCAT.day)]);
            
            % filter for sites with at least 100 measurements
            if length(tem2_lon_lat)>100;
                % filter for sites with temporal spread greater than 10
                % years
                if max((tem2_lon_lat(:,SOCAT.yr)))-min((tem2_lon_lat(:,SOCAT.yr)))>10;
                    if length(unique(tem2_lon_lat(:,SOCAT.yr)))>0.5*(max((tem2_lon_lat(:,SOCAT.yr)))-min((tem2_lon_lat(:,SOCAT.yr))));
                        if length(unique(tem2_lon_lat(:,SOCAT.yr)))>6;
                            k=k+1
                            
                            % remove outliers
                            [B,TF]=rmoutliers(tem2_lon_lat(:,SOCAT.fCO2rec));
                            tem2_lon_lat(TF,:)=[];
                            
%                           % Calculate TA using LIAR (Carter et al., 2018)
                            TA=LIAR(horzcat(tem2_lon_lat(:,[SOCAT.lon SOCAT.lat]),ones(size(tem2_lon_lat,1),1)*3),tem2_lon_lat(:,[SOCAT.SSS SOCAT.SST]),[1,7],'VerboseTF',false); %Add ____,'VerboseTF',false____ if you don't like all of the text on the command line
                            % Add TA_LIAR to SOCAT data
                            tem2_lon_lat(:,SOCAT.TA_Liar)=TA;
                            
                            % Call CO2SYS
                            % **** SYNTAX: [RESULT,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,
                            % **** SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,
                            % **** K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
                            ATADIC=CO2SYS(tem2_lon_lat(:,SOCAT.TA_Liar),tem2_lon_lat(:,SOCAT.fCO2rec),par1type,par2type,tem2_lon_lat(:,SOCAT.SSS),tem2_lon_lat(:,SOCAT.SST),tem2_lon_lat(:,SOCAT.SST),0,0,sil,po4,nh4,h2s,pHscale,k1k2c,kso4c,kfc,tb);
                            
                            % CO2SYS is returning -999 for nan values and
                            % this is screwing everything up AT
                            loc18=find((ATADIC(:,18))<0);
                            ATADIC(loc18,18) = NaN;
                            loc39=find((ATADIC(:,39))<0);
                            ATADIC(loc39,39) = NaN;
                            
                            % From CO2SYS results...
                            % 18 - OmegaAr input
                            % 21 - pH output
                            % 15 - Hfree input
                            % 36 - OmegaAr output
                            % 39 - pH input (Total)
                            % 43 - pH output (Total)
                            tem2_lon_lat(:,SOCAT.OmegaAr)=ATADIC(:,18);
                            tem2_lon_lat(:,SOCAT.pH)=ATADIC(:,39);
                            %tem2(:,SOCAT.Hfree)=ATADIC(:,15);
                            
                            % DETRENDING DATA
                            % Reference year = 2000
                            % tem2(:,SOCAT.fCO2rec_detrended)=tem2(:,SOCAT.fCO2rec)-1.89*(tem2(:,SOCAT.yr)-min(tem2(:,SOCAT.yr)));
                            
                            REF_YEAR = 2000;
                            tem2_lon_lat(:,SOCAT.fCO2rec_detrended)=tem2_lon_lat(:,SOCAT.fCO2rec)-1.89*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
                            tem2_lon_lat(:,SOCAT.pH_detrended)=tem2_lon_lat(:,SOCAT.pH)+0.0018*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
                            tem2_lon_lat(:,SOCAT.OmegaAr_detrended)=tem2_lon_lat(:,SOCAT.OmegaAr)+0.0078*(tem2_lon_lat(:,SOCAT.yr)-REF_YEAR);
                            
                            % Add ordinal date
                            tem2_lon_lat(:,27) = datenum([tem2_lon_lat(:,SOCAT.yr),tem2_lon_lat(:,SOCAT.mon),tem2_lon_lat(:,SOCAT.day)]);
                            
                            % TAKE DAILY AVERAGES -------------------------
                            %
                            % Since this data is meant to be compared to
                            % the monthly GFDL data, and the SOCAT data is
                            % taken every few seconds, we are taking out
                            % some of the shorter-term variability by
                            % taking daily means to make the variability in
                            % this data slightly more comparable to the
                            % monthly GFDL data, without removing all the
                            % variability. 
                            
                            % Step 1. Find all the unique days and prepare
                            % an empty matrix to store the daily means
                            byday = zeros(length(unique(tem2_lon_lat(:,27))),27);
                            % Step 2. Loop through the unique days
                            for p = 1:1:length(unique(tem2_lon_lat(:,12)));
                                A = unique(tem2_lon_lat(:,12));
                                idx = ismember(tem2_lon_lat(:,12),A(p));
                                % Result contains all the data from a
                                % unique day
                                Result = tem2_lon_lat(idx,:);
                                % Step 3. Assign to each row the mean of
                                % the unique days. 
                                % p corresponds to the row
                                % number and unique day. q corresponds to
                                % the number of columns in the original
                                % dataset and in this dataset.
                                for q = 1:1:27;
                                    byday(p,q) = mean(Result(:,q));
                                end
                            end
                            
                            % prep final data 
                            VAR_SOCAT(k,1)=lon;
                            VAR_SOCAT(k,2)=lat;
                            VAR_SOCAT(k,3)=min(tem2_lon_lat(:,SOCAT.yr));
                            VAR_SOCAT(k,4)=max(tem2_lon_lat(:,SOCAT.yr));
                            VAR_SOCAT(k,5)=length(unique(tem2_lon_lat(:,SOCAT.yr)));
%                             VAR_SOCAT(k,6)=nanstd(tem2_lon_lat(:,SOCAT.fCO2rec));
                            VAR_SOCAT(k,6)=nanstd(byday(:,SOCAT.fCO2rec));
%                             VAR_SOCAT(k,7)=nanstd(tem2_lon_lat(:,SOCAT.pH));
                            VAR_SOCAT(k,7)=nanstd(byday(:,SOCAT.pH));
%                             VAR_SOCAT(k,8)=nanstd(tem2_lon_lat(:,SOCAT.OmegaAr));
                            VAR_SOCAT(k,8)=nanstd(byday(:,SOCAT.OmegaAr));
%                             VAR_SOCAT(k,9)=nanstd(tem2_lon_lat(:,SOCAT.fCO2rec_detrended));
                            VAR_SOCAT(k,9)=nanstd(byday(:,SOCAT.fCO2rec_detrended));
%                             VAR_SOCAT(k,10)=nanstd(tem2_lon_lat(:,SOCAT.pH_detrended));
                            VAR_SOCAT(k,10)=nanstd(byday(:,SOCAT.pH_detrended));
%                             VAR_SOCAT(k,11)=nanstd(tem2_lon_lat(:,SOCAT.OmegaAr_detrended));
                            VAR_SOCAT(k,11)=nanstd(byday(:,SOCAT.OmegaAr_detrended));
%                             VAR_SOCAT(k,12) = dist_to_land(k);
%                             VAR_SOCAT(k,12) = mean(tem2_lon_lat(:,27));
%                             VAR_SOCAT(k,13)=min(tem2(:,SOCAT.fCO2rec_detrended));
%                             VAR_SOCAT(k,14)=min(tem2(:,SOCAT.pH_detrended));
%                             VAR_SOCAT(k,15)=min(tem2(:,SOCAT.OmegaAr_detrended));
%                             VAR_SOCAT(k,16)=max(tem2(:,SOCAT.fCO2rec_detrended));
%                             VAR_SOCAT(k,17)=max(tem2(:,SOCAT.pH_detrended));
%                             VAR_SOCAT(k,18)=max(tem2(:,SOCAT.OmegaAr_detrended));
%                     
                        end
                end
                end
            end
    end
end


hearders2={'lon','lat','minYr','maxYr','UniqueYrs','StdfCO2','StdpH','StdOmegaAr','StdfCO2_detrended','StdpH_detrended','StdOmegaAr_detrended'}  ;                 
csvwrite_with_headers('VAR_SOCAT_NanStd.csv',VAR_SOCAT,hearders2);

