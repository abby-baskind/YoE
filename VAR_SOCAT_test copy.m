clear all;
par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par2type =    5; % The first parameter supplied is of type "1", which is "DIC"
sil      =    0; % Concentration of silicate  in the sample (in umol/kg)
po4      =    0; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    13; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; 
load SOCATV2022.mat;
dat=[longitude,latitude,yr,mon,day,hh,mm,ss,SST,sal,fCO2rec];
headers={'lon','lat', 'yr' 'mon','day','hh','mm','ss','SST','SSS','fCO2rec','TA_Liar','OmegaAr','OmegaCa','pH','Hfree','fCO2rec_detrended','OmegaAr_detrended','pH_detrended'};
 for n=1:size(headers,2)
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
                if max((tem2(:,SOCAT.yr)))-min((tem2(:,SOCAT.yr)))>10
                    if length(unique(tem2(:,SOCAT.yr)))>0.5*(max((tem2(:,SOCAT.yr)))-min((tem2(:,SOCAT.yr))))
                        if length(unique(tem2(:,SOCAT.yr)))>6
                            k=k+1;
                            [B,TF]=rmoutliers(tem2(:,SOCAT.fCO2rec));
                            tem2(TF,:)=[];
                            TA=LIAR(horzcat(tem2(:,[SOCAT.lon SOCAT.lat]),ones(size(tem2,1),1)*3),tem2(:,[SOCAT.SSS SOCAT.SST]),[1,7],'VerboseTF',false); %Add ____,'VerboseTF',false____ if you don't like all of the text on the command line 
                            tem2(:,SOCAT.TA_Liar)=TA;
                            ATADIC=CO2SYS(tem2(:,SOCAT.TA_Liar),tem2(:,SOCAT.fCO2rec),par1type,par2type,tem2(:,SOCAT.SSS),tem2(:,SOCAT.SST),tem2(:,SOCAT.SST),0,0,sil,po4,pHscale,k1k2c,kso4c);
                            tem2(:,SOCAT.OmegaAr)=ATADIC(:,31);
                            tem2(:,SOCAT.pH)=ATADIC(:,3);
                            %tem2(:,SOCAT.Hfree)=ATADIC(:,13);
                            % chanbe from min year to a reference year
                            tem2(:,SOCAT.fCO2rec_detrended)=tem2(:,SOCAT.fCO2rec)-1.89*(tem2(:,SOCAT.yr)-min(tem2(:,SOCAT.yr)));
                            tem2(:,SOCAT.pH_detrended)=tem2(:,SOCAT.pH)+0.0018*(tem2(:,SOCAT.yr)-min(tem2(:,SOCAT.yr)));
                            tem2(:,SOCAT.OmegaAr_detrended)=tem2(:,SOCAT.OmegaAr)+0.0078*(tem2(:,SOCAT.yr)-min(tem2(:,SOCAT.yr)));
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
                    
                        end
                end
                end
            end
    end
end

hearders2={'lon','lat','minYr','maxYr','UniqueYrs','StdfCO2','StdpH','StdOmegaAr','StdfCO2_detrended','StdpH_detrended','StdOmegaAr_detrended'}  ;                 
csvwrite_with_headers('VAR_SOCAT_NanStd.csv',VAR_SOCAT,hearders2);
 