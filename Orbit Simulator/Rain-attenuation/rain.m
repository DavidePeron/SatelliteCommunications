% Estimation of rain attenuation following the simplified model

lat = input('Insert the latitude: ');
lon = input('Insert the longitude:');

if lon < 0
    lon = 360 + lon;
end

%elevation=input('Insert the elevation (theta): ')
% lat=61.267865; %latitude and longitude 
% lon=96.608223;
elevation= 43.168; %best elevation 
%pol=input('Insert the polarization V or H: ,'s')

p=0.01; %rainfall rate exceeded for p% of an average year
Re=8500;
[RainRate,p0] = Rec837_5(p,lat,lon);% RainRate exceeded for 0.01% of the time for the ground station segment

%hr: Height of the rain for the 0.01% of the time
if(lat<20)
    rophi=0.6; %rophi is a empiric factor of height reduction
elseif (lat>=20)&&(lat<40)
    rophi=0.6+0.02*(lat-20);
else
    rophi=1;
end
hr=rophi*[5.1-2.15*log10(1+10^((lat-27)/25))]; %km

%hs: Height of the GS above sea level REC ITU-R P.1511
load TOPOLAT.txt
load TOPOLON.txt
load TOPO_0DOT5.txt

temp_lat = lat - floor(abs(lat))*sign(lat);
temp_lon = lon - floor(abs(lon))*sign(lon);

if abs(temp_lat) < 0.5
    lat = lat - temp_lat;
elseif (abs(temp_lat) >= 0.5)&&(abs(temp_lat) < 0.75)
    lat = lat - temp_lat + 0.5*sign(lat);
elseif abs(temp_lat) >= 0.75
    lat = lat - temp_lat + 1*sign(lat);
end

if abs(temp_lon) < 0.5
    lon = lon - temp_lon;
elseif (abs(temp_lon) >= 0.5)&&(abs(temp_lon) < 0.75)
    lon = lon - temp_lon + 0.5*sign(lon);
elseif abs(temp_lon) >= 0.75
    lon = lon - temp_lon + 1*sign(lon);
end

%Now we need to round the latitude and the longitude to the first decimal
%and still keeping the variables as double, otherwise the find function
%will never find a match with the values of the txt file, which have to be
%double

% % % temp_lat = lat - floor(abs(lat))*sign(lat);
% % % temp_lon = lon - floor(abs(lon))*sign(lon);
% % % 
% % % if abs(temp_lat) < 0.5
% % %     lat = lat - temp_lat;
% % % elseif (abs(temp_lat) >= 0.5)&&(abs(temp_lat) < 0.75)
% % %     lat = lat - temp_lat + 0.5*sign(lat);
% % % elseif abs(temp_lat) >= 0.75
% % %     lat = lat - temp_lat + 1*sign(lat);
% % % end
% % % 
% % % if abs(temp_lon) < 0.5
% % %     lon = lon - temp_lon;
% % % elseif (abs(temp_lon) >= 0.5)&&(abs(temp_lon) < 0.75)
% % %     lon = lon - temp_lon + 0.5*sign(lon);
% % % elseif abs(temp_lon) >= 0.75
% % %     lon = lon - temp_lon + 1*sign(lon);
% % % end

[x,y] = find((lat==TOPOLAT)&(lon==TOPOLON));

hs=TOPO_0DOT5(x,y); %km

%Ls: Length of the path throught the rain

if (elevation >5)
    Ls=((hr-hs)/(sin(degtorad(elevation))));
else
    Ls=(2*(hr-hs))/((sqrt(sin(degtorad(elevation)))^2+(2*(hr-hs))/Re)+sin(degtorad(elevation)));
end

%rp: factor of no homogeneity of the rain
rp = 90/(90+4*Ls*cos(degtorad(elevation)));

% ITU-R P.838-3
fUL=14.414; %Ghz
fDL=11.114; %GHz
kUL=(4.21e-5)*(fUL^2.42);
kDL=(4.21e-5)*(fDL^2.42);
alphaUL=(1.41)*(fUL^-0.0779);
alphaDL=(1.41)*(fDL^-0.0779);

%Total rain attenuation exceeded the 0.01% of the time 
AtotUL=kUL*(RainRate^alphaUL)*Ls*rp; %dB
AtotDL=kDL*(RainRate^alphaDL)*Ls*rp;





