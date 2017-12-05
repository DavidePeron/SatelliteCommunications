function [RainRate,p0] = Rec837_5(p,lat,lon)
%
%-------------------------------------------------------------------------------------------------------------------
%
% function [RainRate,p0] = Rec837_5(p,lat,lon)
%
% This function calculates the rainfall rate exceeded for p% of an average year, 
% according to the prediction method defined in Recommendation ITU-R P.837-5 (Annex 1).
%
% Input parameters :
%   p   = probability level (%)     (may be a vector)
%   lat = latitude          (deg)   (may be a vector) 
%   lon = longitude         (deg)   (may be a vector)
%
% Output parameters :   
%   rr = rainfall rate exceeded for p% of an average year  (mm/h)
%   p0 = percentage probability of rain in an average year (%)
%
% Note :
%   The function automatically loads meteorological input parameters required by 
%   Rec. ITU-R P.837-5 from the following files :
%       ESARAIN_MT_v5.TXT :     mean annual rainfall amount                   (mm)
%       ESARAIN_BETA_v5.TXT :   ratio of convective to total rainfall amount  (-)
%       ESARAIN_PR6_v5.TXT :    probability of rainy 6-hours periods          (%)
%       ESARAIN_LAT_v5.TXT :    latitudes of grid points (90:-1.125:-90)      (deg)
%       ESARAIN_LON_v5.TXT :    longitudes of grid points (0:1.125:360)       (deg)
%
% Called functions :
%   No
%
% Necessary toolboxes :
%   No
%
% Kevin Paulson, University of Hull
% based on a programme by Giulio BLARZINO,
% ONERA, France
% Any questions : email <k.paulson@hull.ac.uk> 
%
% rel. 0.0
% release history:
%   0.0 (1/4/2011)    - original version
%
%-------------------------------------------------------------------------------------------------------------------
%
%           **************************
%           ** Check Input Data     **
%           **************************

%   Check latitude and longitude arrays are the same size
if sum( size(lon) ~= size(lon) ) ~=0
    error('The lat and lon vectors must be the same size');
end

%   Ensure all longitudes are between 0 and 360
lon = mod(lon,360);

%
%           **********************************
%           ** Input Climate Parameters     **
%           **********************************

% load input meteorological paramaters
load ESARAIN_LAT_v5.TXT -ascii ; lat_e40 = ESARAIN_LAT_v5 ;
load ESARAIN_LON_v5.TXT -ascii ; lon_e40 = ESARAIN_LON_v5 ;
load ESARAIN_MT_v5.TXT -ascii ; mt = ESARAIN_MT_v5 ;
load ESARAIN_BETA_v5.TXT -ascii ; conv_ratio = ESARAIN_BETA_v5 ;
load ESARAIN_PR6_v5.TXT -ascii ; pr6 = ESARAIN_PR6_v5 ;

%
%           *****************************
%           ** Perform calculation     **
%           *****************************


% bi-linear interpolation of parameters @ the required coordinates
pr6i  = interp2(lon_e40,lat_e40,pr6,lon,lat,'linear');
mti   = interp2(lon_e40,lat_e40,mt,lon,lat,'linear');
betai = interp2(lon_e40,lat_e40,conv_ratio,lon,lat,'linear');

% extract mean annual rainfall amount of stratiform type
msi = mti.*(1 - betai);

% percentage probability of rain in an average year
p0 = pr6i.*(1 - exp(-0.0079.*(msi./pr6i)));

% Loop over each (lat,lon) point and calculate rain rate exceedance

nPoint = numel(lat);            % number of (lat,lon) points
nProb = numel(p);               % number of exceedances
RainRate = zeros(nProb,nPoint); % Declare array to hold rain rates

for iPoint = 1:nPoint
    
    rr = zeros(numel(p),1);     % Array to hold rain rates for this point
    if isnan(p0(iPoint))        % catch the case where P0 is not defined
       p0(iPoint) = 0;
    else
        % P0 is defined so do the calculation
        ThisP0 = p0(iPoint);
       ix = find(p > ThisP0);
       if ~isempty(ix),
          rr(ix) = 0;
       end;
       ix = find(p <= ThisP0);
       if ~isempty(ix),
          a = 1.09;
          b = mti(iPoint)/(21797*ThisP0);
          c = 26.02*b;
          A = a*b;
          B = a + c*log(p(ix)/ThisP0);
          C = log(p(ix)/ThisP0);
          rr(ix) = (-B + sqrt(B.^2 - 4*A*C))/(2*A);
       end
    end
    RainRate(:,iPoint) = rr;

end
