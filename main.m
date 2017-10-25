clear all;
close all;
%% Constant parameters
mu = 3.986004418e14; %Kepler constant [m^3/s^2]
re = 6377000; %Earth radius [m] 
Go = 1.727564365843028; % (rad) Distance in rad of Greenwich meridian from vernal equinox (X-axis)
n_rev =2; %Number of revolutions
runspeed = 300;
omega_earth = 7.292115855377074e-005; % (rad/sec) 
El = 15; %Elevation in degrees
%% Some parameters about orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TUNDRA %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
inc = deg2rad(63.4); %inclination
T = 86400; %period in sec (86400s = 24h)
a = 42164000; % altitude [m]
e = 0.25; %Eccentricity
w = deg2rad(270); %Argument of perigee [deg]
raan = deg2rad(280); %Right ascension of the ascending node [deg]

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOLNIYA %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% inc = deg2rad(63.4); %inclination
% T = 43200; %period in sec (86400s = 24h)
% a = 26556000; % altitude [m]
% e = 0.71; %Eccentricity
% w = deg2rad(270); %Argument of perigee [deg]
% raan = deg2rad(30); %Right ascension of the ascending node [deg]

%WGS84 ellipsoid constants
a_WGS84 = 6378137;
e_WGS84 = 8.1819190842622e-2;
b = sqrt(a_WGS84^2*(1-e_WGS84^2)); %Semi-minor axis
ep = sqrt((a_WGS84^2 - b^2)/b^2); %

eta_0 = 2*pi/T; %Angular velocity of the fictitious satellite (rad/sec)

%Initial position (focal distance)
phi = 0:0.01:2*pi;
r = zeros(1,length(phi));
r(1,:) = a*(1-e^2)./(1+e*cos(phi(1,:))); % r(phi) position of the satellite along the orbit

t = 0:runspeed:T*n_rev; % time [s]
t_p = 0:T:T*n_rev; % time at perigee pass at each cycle

%Transformation Matrix from orbital to ECI Coordinates

orbital_to_ECI = [cos(raan)*cos(w) - sin(raan)*cos(inc)*sin(w), -cos(raan)*sin(w) - sin(raan)*cos(inc)*cos(w), sin(raan)*sin(inc);
     sin(raan)*cos(w) + cos(raan)*cos(inc)*sin(w), -sin(raan)*sin(w) + cos(raan)*cos(inc)*cos(w), -cos(raan)*sin(inc);
     sin(inc)*sin(w), sin(inc)*cos(w), cos(inc)];

% Compute Mean Anomaly and Eccentric Anomaly
M = zeros(1,length(t));
E = zeros(1,length(t));
r_t = zeros(1,length(t));
v_t = zeros(1,length(t));
x_0 = zeros(1,length(t));
y_0 = zeros(1,length(t));
z_0 = zeros(1,length(t));
phi_t = zeros(1,length(t));
ECI_coord = zeros(3,length(t));
ECI_velocity = zeros(4,length(t));
ECEF_coord = zeros(3,length(t));
long = zeros(1,length(t));
lat = zeros(1,length(t));

%% Earth-satellite geometry parameters
elong_ext = zeros(1,length(t)); %longitude of the extreme toward east
nlat_ext = zeros(1,length(t)); %latitude of the extreme toward north
wlong_ext = zeros(1,length(t)); %longitude of the extreme toward west
slat_ext = zeros(1,length(t)); %latitude of the extreme toward south
lambda = zeros(1,length(t)); %vector for all the lambda values
rho = zeros(1,length(t));
lambda_0 = zeros(1,length(t));
nadir = zeros(1,length(t));

%% Plot of the static components

%Initializing the Drawing Space
set(gcf,'Menubar','default','Name','Orbit Visualization', ... 
    'NumberTitle','off','Position',[10,350,750,750]); 
lim=(1+e)*a;%Setting the limits of the graph
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
title('Orbital Visualization');

%Plotting the Earth
%Julian Date of 10 January 2017 at noon
A = 2017; %Year
DTA = 10; %Days from 1 Jan of year A
NAB1900 = 29; %Number of leap years since 1900
TU = 12; %Hours in the day
JD = 2415020 + 365*(A - 1900) + DTA + NAB1900 + TU/24 - 0.5;

T_c = (JD - 2415020)/36525;
alpha_go = 99.6909833 + 36000.7689*T_c + 3.8708e-4*T_c^2;

%Number of Days since Jan 1, 2000
% J2000_days=155727/24; % = 4th October, 2017 http://www.timeanddate.com/counters/year2000.html
equat_rad=6378137.00;
polar_rad=6356752.3142;
[xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
load('topo.mat','topo','topomap1');
topo2 = [topo(:,181:360) topo(:,1:180)];
pro.FaceColor= 'texture';
pro.EdgeColor = 'none';
pro.FaceLighting = 'phong';
pro.Cdata = topo2;
earth= surface(xx,yy,zz,pro);
colormap(topomap1);

rotate (earth, [0 0 1], 0);
Xaxis= line([-1e7 1e7],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis= line([0 0],[-1e7 1e7],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
rotate (Xaxis, [0 0 1], 0);
rotate (Yaxis, [0 0 1], 0);
Sun=light('Position',[1 0 0],'Style','infinite');

%Plotting the ECI Axes
line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')

% %Plotting the initial poisition of the satellite
% plot3 (ECI_coord(1,1), ECI_coord(2,1), ECI_coord(3,1),'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','green','MarkerSize', 10);
% 
% %Plotting Initial Velocity Vector
% line([ECI_coord(1,1) ECI_coord(1,1)+ECI_velocity(1,1)],[ECI_coord(2,1) ECI_coord(2,1)+ECI_velocity(3,1)],[ECI_coord(3,1) ECI_coord(3,1)+ECI_velocity(4,1)],'Color', 'red','Marker','.','LineWidth', 1.5, 'MarkerSize', 8,'LineStyle','-');

%% DRAWING THE DYNAMIC VISUALIZATION COMPONENTS%

%  sphere_position = cell(1,length(t)); %1,length(t)
%  position = cell(1,length(t)); %1,length(t)

for i = 1:length(t)

    flag = 0;
    %Mean Anomaly
    M(1,i) = eta_0*(t(i) - t_p(1,floor(t(i)/T)+1));
    
    % Computation of Eccentric Anomaly with Newton-Raphson Method
    E(i) = pi;
    
    %Newton-Raphson Method
    while flag == 0
        temp_1 = E(i);
        E(i) = E(i) + (M(i) + e*sin(E(i)) - E(i))/(1 - e*cos(E(i)));
        temp_2 = E(i);
        if abs(temp_1 - temp_2) < 1e-2
            flag = 1;
        end
    end
    
    
    r_t(i) = a*(1-e*cos(E(i)));
    v_t(i) = sqrt(mu*(2/r_t(i)-1/a));
    phi_t(i) = acos(1/e*((a*(1-e^2))/r_t(i)-1));
    
    if(E(i) > pi)
        phi_t(i) = 2*pi - phi_t(i);
    end
    
    x_0(i) = r_t(i)* cos(phi_t(i));
    y_0(i) = r_t(i)* sin(phi_t(i));
    ECI_coord(:,i) = orbital_to_ECI*[x_0(i); y_0(i); z_0(i)];
    
    rotate (earth, [0 0 1], (360/length(t)), [0,0,0]);
    rotate (Xaxis, [0 0 1], (360/length(t)), [0,0,0]);
    rotate (Yaxis, [0 0 1], (360/length(t)), [0,0,0]);

    %Drawing the red sphere
    sphere_position(i)=plot3 (ECI_coord(1,i), ECI_coord(2,i), ECI_coord(3,i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
    position(i)=line([0 ECI_coord(1,i)],[0 ECI_coord(2,i)], [0 ECI_coord(3,i)],'Color', 'yellow', 'LineWidth', 2);
    
    if(i~=1)
        set (sphere_position(i-1), 'Visible', 'off');
        set (position(i-1), 'Visible', 'off');
    end
    if (i~=1 && i<=length(t))
        line([ECI_coord(1,i-1) ECI_coord(1,i)],[ECI_coord(2,i-1) ECI_coord(2,i)], [ECI_coord(3,i-1) ECI_coord(3,i)],'Color', 'black', 'LineWidth', 1.5);
    end
    
    %Compute ECEF Coordinates
    %Transformation Matrix from ECI to ECEF Coordinates

    theta = alpha_go + 0.25068447*t(i)/60;
    theta = wrapTo2Pi(deg2rad(theta));
    ECI_to_ECEF = [cos(theta) sin(theta) 0; 
                  -sin(theta) cos(theta) 0; 
                       0          0      1];
    ECEF_coord(:,i) = ECI_to_ECEF*ECI_coord(:,i);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transformation to LLA with formulas from slides does not work !!!
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = sqrt(ECEF_coord(1,i).^2 + ECEF_coord(2,i).^2);
    th = atan2(a_WGS84*ECEF_coord(3,i), b*p); %theta
    long(i) = atan2(ECEF_coord(2,i), ECEF_coord(1,i));
    lat(i) = atan2(ECEF_coord(3,i)+ep^2.*b.*sin(th).^3 , p-e_WGS84^2.*a_WGS84.*cos(th).^3);
    
    N = a_WGS84./sqrt(1-e_WGS84.^2*sin(lat(i)).^2);
    alt = p./(cos(lat(i))-N);
    
    long(i) = rad2deg(long(i));
    lat(i) = rad2deg(lat(i));

    %return long in [0,360] range
    long(i) = mod(long(i),360);
    
    %Pause
    pause (0.01);
    
end

%% Ground Track

perigee=a*(1-e);
apogee=a*(1+e);
altitude_low=(perigee-re)/1000;
altitude_high=(apogee-re)/1000;

if altitude_low>200 && altitude_high<2000
    'Low Earth Orbit'
end
if altitude_low>2000 && altitude_high<35786
    'Middle Earth Orbit'
end

figure (2);
set(gcf,'Menubar','none','Name','Earth Track', ... 
    'NumberTitle','off','Position',[10,350,1000,500]); 
hold on
image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
colormap(topomap1);
axis equal
axis ([0 360 -90 90]);
h1 = plot(0,0,'o','Markersize',1);

%Ground Stations plot
% plot (167.717,8.717,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
% plot (360-76.496, 42.440,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
for i=1:length(t)
    plot (long(i),lat(i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 2);
    if (i~=1 && abs(long(i-1)-long(i))<100)
        line([long(i-1) long(i)],[lat(i-1) lat(i)],'Color', 'red', 'LineWidth', 2);
    end
    
    %% Nadir angle estimation
    rho(i) = asin(re/r_t(i));
    lambda_0(i) = pi/2 - rho(i);
    beta = asin(cos(deg2rad(El))*sin(rho(i)));
    nadir(i) = beta;

    %% Lambda and arch estimation 
    lambda(i) = 90 - rad2deg(nadir(i)) - El; %central angle in degrees
    
    arch = (2*pi*r_t(i)/360)*lambda(i); %arch between subsat. point and extreme position of a user
    l = 2*re*sin(deg2rad(lambda(i)/2)); %segment between subsat. point and extreme position of a user
    
    %Finding the LLA coordinates of the extreme position of a user
    %We consider positive angles towards east and north and negative ones 
    %towards west and south; 
    elong_ext(i) = (long(i) + lambda(i)) - 360*floor((long(i) + lambda(i))/360);
    wlong_ext(i) = (long(i) - lambda(i));
    nlat_ext(i) = (lat(i) + lambda(i)); %- 90*floor(lat(i)+lambda(i)/90);
    slat_ext(i) = lat(i) - lambda(i);
    
    delete(h1);
    h1 = plot(long(i),lat(i),'o','Markersize',lambda(i)*2);
    
    pause (0.001);
end

%% 
% %Plot of r(phi)
% figure(3);
% set(gcf,'Menubar','default','Name','r(phi)'); 
% plot(phi,r);
% 
% %Plot of r(t)
% figure(4);
% set(gcf,'Menubar','default','Name','r(t)'); 
% plot(t, r_t);
% 
%Plot of latitude
figure(5);
set(gcf,'Menubar','default','Name','LLA Coordinates and Parameters'); 
subplot(2,2,1);
plot(t, lat);
title('Latitude(t)','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label

%Plot of longitude
subplot(2,2,2);
plot(t, long);
title('Longitude(t)','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label
% 
% %Plot of Anomalies
% figure(7);
% set(gcf,'Menubar','default','Name','Anomalies'); 
% plot(t, M, 'r');
% hold on
% plot(t, E, 'b');
% hold on
% plot(t, phi_t, 'g');

%Plot of lambda (radius of the circle)
subplot(2,2,3);
plot(t, lambda);
title('$$Central ~angle ~\lambda(t)$$','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label

%Time of visibility, max Angular Speed and Azimuth range [s]
Tv=(T/pi)*acos(cos(deg2rad(max(lambda)))/cos(deg2rad(min(lambda))));
AngularSpeed = zeros(1,length(t));
Dmin=re*(sin(deg2rad(min(lambda)))/sin(min(nadir)))
AngularSpeed(1,:)=(2*pi*(r_t(1,:)))/(T*Dmin)
maxAngularSpeed=max(AngularSpeed)
AzimuthRange=2*acos(tan(min(deg2rad(lambda)))/(tan(max(deg2rad(lambda)))))
%%% Number of satellites 
gamma=deg2rad(max(lambda)); %central angle
s=2*pi/(sqrt(3)*gamma); %number of satellites per orbiltal plane
s=floor(s);
orbitalp=2*pi/(3*gamma); %number of orbital planes
orbitalp=ceil(orbitalp);
totalnsatellites=s*orbitalp; %total number of satellites
sepsat=360/s;%separation of satellites in each orbital plane
relativephasing=sepsat/2; %relative phasing between satellites in adjacent planes
GAMMA=rad2deg(acos(cos(gamma)/(cos(pi/s))));%hald of ground swath width
alpha=GAMMA+rad2deg(gamma)%separation between orbital planes

lat_deg = lat;
long_deg = long;
%% Azimuth estimation
prompt = '\n Insert a latitude for the ground station: ';
lat_gs = deg2rad(input(prompt));
prompt = '\n Insert a longitude between [-180; +180] for the ground station: ';
long_gs = deg2rad(input(prompt));
azimuth = zeros(1,length(t));
elevation = zeros(1,length(t));
% Save latitude and longitude in degree to future checks

for i = 1:length(t)
    if long(i) > 180
        long(i) = long(i) - 360;
    end
    long(i) = deg2rad(long(i));
    lat(i) = deg2rad(lat(i));
%     L = long(i) - long_gs;
%     cosine_phi = cos(L)*cos(lat_gs);
%     a_az = asin(sin(L)/sin(cosine_phi));
%     if lat(i) > 0 && L > 0
%         azimuth(i) = 180 - rad2deg(a_az);
%     end
%     if lat(i) > 0 && L < 0
%         azimuth(i) = 180 + rad2deg(a_az);
%     end
%     if lat(i) < 0 && L > 0
%         azimuth(i) = rad2deg(a_az);
%     end
%     if lat(i) < 0 && L < 0
%         azimuth(i) = 360 - rad2deg(a_az);
%     end      
    cos_lambda = cos(lat_gs)*cos(lat(i))*cos(long_gs - long(i)) + sin(lat_gs)*sin(lat(i));
    elevation(i) = rad2deg(acos(sin(acos(cos_lambda))/sqrt(1+(re/r_t(i))^2+2*re/r_t(i)*cos_lambda)));
end
% 
% % %Plot of Azimuth 
% % figure(9);
% % set(gcf,'Menubar','default','Name','azimuth(t)'); 
% % plot(t, azimuth);
% % title('$$Azimuth(t)$$','interpreter','latex');
% % xlabel('Time'); % x-axis label
% % ylabel('Degree'); % y-axis label
% 
%Plot of Elevation 
subplot(2,2,4); 
plot(t, elevation);
title('$$Elevation(t)$$','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label

pks = findpeaks(long_deg(1:ceil(length(long_deg)/2))); %pks = most eastern longitude
max_long = min(pks); %the smaller extreme in longitude
locs = find(long_deg == max_long); %index of max_long in the long vector

ext_lat = lat_deg(locs); %ext_lat = equivalent latitude at max_long
lat_index = find(lat_deg == ext_lat); %lat_index = index of the most western longitude
if lat_index(1) == locs
    lat_locs = lat_index(2);
else
    lat_locs = lat_index(1);
end
min_long = long_deg(lat_locs); %most western longitude
radius = lambda(locs); %radius of the circle of coverage for the max_long point
rext_east = max_long + radius;
lext_east = min_long + radius;
semi_cov = abs(rext_east - lext_east);

cov = 2*radius; %actual double coverage for tundra orbit in the north part
n_orbits = ceil(200/cov);

