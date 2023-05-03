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
load SOCATV2022.mat;
dat=[longitude,latitude,yr,mon,day,hh,mm,ss,SST,sal,fCO2rec];
headers={'lon','lat', 'yr' 'mon','day','hh','mm','ss','SST','SSS','fCO2rec','TA_Liar','OmegaAr','OmegaCa','pH','Hfree','fCO2rec_detrended','OmegaAr_detrended','pH_detrended', 'dist_to_land','fCO2rec_detrended_min','OmegaAr_detrended_min','pH_detrended_min','fCO2rec_detrended_max','OmegaAr_detrended_max','pH_detrended_max'};
 for n=1:size(headers,2);
    SOCAT.(headers{1,n})=n;
 end  
 VAR_SOCAT=[];
 k=0
for lon=0:1:360;
        loc1=find((dat(:,SOCAT.lon))<lon+1 & (dat(:,SOCAT.lon))>=lon);
        tem1=dat(loc1,:);
        for lat=-90:1:90;
            loc2=find((tem1(:,SOCAT.lat))<lat+1 & (tem1(:,SOCAT.lat))>=lat);
            tem2=tem1(loc2,:);
            if length(tem2)>100;
                if max((tem2(:,SOCAT.yr)))-min((tem2(:,SOCAT.yr)))>10;
                    if length(unique(tem2(:,SOCAT.yr)))>0.5*(max((tem2(:,SOCAT.yr)))-min((tem2(:,SOCAT.yr))));
                        if length(unique(tem2(:,SOCAT.yr)))>6;
                            k=k+1
                            [B,TF]=rmoutliers(tem2(:,SOCAT.fCO2rec));
                            tem2(TF,:)=[];
                            
%                           % Calculate TA using LIAR (Carter et al., 2018)
                            TA=LIAR(horzcat(tem2(:,[SOCAT.lon SOCAT.lat]),ones(size(tem2,1),1)*3),tem2(:,[SOCAT.SSS SOCAT.SST]),[1,7],'VerboseTF',false); %Add ____,'VerboseTF',false____ if you don't like all of the text on the command line
                            % Add TA_LIAR to SOCAT data
                            tem2(:,SOCAT.TA_Liar)=TA;
                            
                            % Call CO2SYS
                            % **** SYNTAX: [RESULT,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,
                            % **** SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,
                            % **** K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
                            ATADIC=CO2SYS(tem2(:,SOCAT.TA_Liar),tem2(:,SOCAT.fCO2rec),par1type,par2type,tem2(:,SOCAT.SSS),tem2(:,SOCAT.SST),tem2(:,SOCAT.SST),0,0,sil,po4,nh4,h2s,pHscale,k1k2c,kso4c,kfc,tb);
                            
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
                            tem2(:,SOCAT.OmegaAr)=ATADIC(:,18);
                            tem2(:,SOCAT.pH)=ATADIC(:,39);
                            %tem2(:,SOCAT.Hfree)=ATADIC(:,15);
                            
                            % DETRENDING DATA
                            % Reference year = 1980
                            % tem2(:,SOCAT.fCO2rec_detrended)=tem2(:,SOCAT.fCO2rec)-1.89*(tem2(:,SOCAT.yr)-min(tem2(:,SOCAT.yr)));
                            tem2(:,SOCAT.fCO2rec_detrended)=tem2(:,SOCAT.fCO2rec)-1.89*(tem2(:,SOCAT.yr)-1980);
                            tem2(:,SOCAT.pH_detrended)=tem2(:,SOCAT.pH)+0.0018*(tem2(:,SOCAT.yr)-1980);
                            tem2(:,SOCAT.OmegaAr_detrended)=tem2(:,SOCAT.OmegaAr)+0.0078*(tem2(:,SOCAT.yr)-1980);
                            VAR_SOCAT(k,1)=lon;
                            VAR_SOCAT(k,2)=lat;
                            VAR_SOCAT(k,3)=min(tem2(:,SOCAT.yr));
                            VAR_SOCAT(k,4)=max(tem2(:,SOCAT.yr));
                            VAR_SOCAT(k,5)=length(unique(tem2(:,SOCAT.yr)));
                            VAR_SOCAT(k,6)=nanstd(tem2(:,SOCAT.fCO2rec));
                            VAR_SOCAT(k,7)=nanstd(tem2(:,SOCAT.pH));
                            VAR_SOCAT(k,8)=nanstd(tem2(:,SOCAT.OmegaAr));
                            VAR_SOCAT(k,9)=nanstd(tem2(:,SOCAT.fCO2rec_detrended));
                            VAR_SOCAT(k,10)=nanstd(tem2(:,SOCAT.pH_detrended));
                            VAR_SOCAT(k,11)=nanstd(tem2(:,SOCAT.OmegaAr_detrended));
                            VAR_SOCAT(k,12) = dist_to_land(k);
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


hearders2={'lon','lat','minYr','maxYr','UniqueYrs','StdfCO2','StdpH','StdOmegaAr','StdfCO2_detrended','StdpH_detrended','StdOmegaAr_detrended', 'dist_to_land', 'fCO2rec_detrended_min','OmegaAr_detrended_min','pH_detrended_min','fCO2rec_detrended_max','OmegaAr_detrended_max','pH_detrended_max'}  ;                 
csvwrite_with_headers('VAR_SOCAT_NanStd.csv',VAR_SOCAT,hearders2);

