clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EARTH'S CONSTANTS USED IN THIS SIMULATOR
% mu: Kepler constant [m^3/s^2]
% re: Earth radius [m]
% Go: Distance in of Greenwich meridian from vernal equinox (X-axis) [rad]
% omega_earth: Angular velocity of the earth [rad/sec]
% a_WGS84: Major semi-axis of the earth modelled as the ellipsoid WGS84
% e_WGS84: Parameter of the earth modelled as the ellipsoid WGS84
% b: Minor semi-axis of the earth modelled as the ellipsoid WGS84
% ep: Parameter of the earth modelled as the ellipsoid WGS84
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('earth_constants.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORBIT'S CONSTANTS
% inc: inclination of the orbit [rad]
% T: orbital period [s]
% a: altitude measured from the center of the earth [m]
% e: orbit's eccentricity
% w: argument of the perigee [rad]
% raan: Right Ascension of the ascending node [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('tundra.mat');

%%%%%% To create a personalized orbit, decomment the following code %%%%%%%
% inc = deg2rad(63.4); %inclination
% T = 43200; %period in sec (86400s = 24h)
% a = 26556000; % altitude [m]
% e = 0.71; %Eccentricity
% w = deg2rad(270); %Argument of perigee [deg]
% raan = deg2rad(30); %Right ascension of the ascending node [deg]

% Simulation parameters
n_rev =2; %Number of revolutions
runspeed = 300; % Speed of the simulation
El = 15; %Minimum elevation in degrees
n_sat = 3;
color = ['r','g','c','y','w'];
raan = deg2rad(25);
lat_gs = [80, 80];
long_gs = [-180, 0];
% raan_sat = [raan, raan + deg2rad(-242.4783429000001), raan + deg2rad(-122.1497972999974)];
raan_sat = zeros(1,n_sat)+raan; % Initial raan vector
t = 0:runspeed:T*n_rev; % time [s]

r = cell(n_sat); %cell vector containing the coordinates of the three satellites
phi = cell(n_sat);

%% RAAN ESTIMATION FOR DIFFERENT SATELLITE IN THE SAME ORBITAL PLANE%

%Computation of polar coordination
[r_t, phi_t] = get_polar(t, T, e, a, mu, n_rev);

for i=0:n_sat-1    
    r{i+1} = circshift(r_t,(T/n_sat*i)./runspeed,2);
    phi{i+1} = circshift(phi_t,(T/n_sat*i)./runspeed,2);
%     raan_sat(i+1) = raan_sat(i+1) + phi{i+1}(1);
end

for i=2:n_sat
    raan_sat(i) = raan_sat(i) + deg2rad(360/n_sat*(i-1));
end

% [~, long, ~] = polar_to_LLA(t, r, phi, inc, w, raan_sat);
% 
% for i=2:n_sat
%     long_temp = circshift(long{1},(T/n_sat*(i-1))./runspeed,2);
%     raan_diff = long{1}(1) - long{i}(1);
%     raan_sat(i) = raan_sat(i) + deg2rad(raan_diff);
% end

%% COMPUTATION OF THE TRAJECTORY%

[ECI, long, lat] = polar_to_LLA(t, r, phi, inc, w, raan_sat);

%% 3-D ANIMATION

%Initializing the Drawing Space and static components
earthmap = imread('planisphere2.jpg');

earthplot3D(earthmap, t, r, re, El, e, a, ECI, color);

%% Ground Track

groundtrack_plot(earthmap, t, r, re, El, long, lat, color, lat_gs, long_gs);

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
plot(t, lat{1});
title('Latitude(t)','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label
grid on;

%Plot of longitude
subplot(2,2,2);
plot(t, long{1});
title('Longitude(t)','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label
grid on;

% % % %Time of visibility, max Angular Speed and Azimuth range [s]
% % % Tv=(T/pi)*acos(cos(deg2rad(max(lambda)))/cos(deg2rad(min(lambda))));
% % % AngularSpeed = zeros(1,length(t));
% % % Dmin=re*(sin(deg2rad(min(lambda)))/sin(min(nadir)))
% % % AngularSpeed(1,:)=(2*pi*(r_t(1,:)))/(T*Dmin)
% % % maxAngularSpeed=max(AngularSpeed)
% % % AzimuthRange=2*acos(tan(min(deg2rad(lambda)))/(tan(max(deg2rad(lambda)))))
% % % %%% Number of satellites 
% % % gamma=deg2rad(max(lambda)); %central angle
% % % s=2*pi/(sqrt(3)*gamma); %number of satellites per orbiltal plane
% % % s=floor(s);
% % % orbitalp=2*pi/(3*gamma); %number of orbital planes
% % % orbitalp=ceil(orbitalp);
% % % totalnsatellites=s*orbitalp; %total number of satellites
% % % sepsat=360/s;%separation of satellites in each orbital plane
% % % relativephasing=sepsat/2; %relative phasing between satellites in adjacent planes
% % % GAMMA=rad2deg(acos(cos(gamma)/(cos(pi/s))));%hald of ground swath width
% % % alpha=GAMMA+rad2deg(gamma)%separation between orbital planes

%% Azimuth estimation

for k = 1:length(lat_gs)

    [azimuth, elevation, best_elevation, best_azimuth] = view_angles(t, lat_gs(k), long_gs(k), long, lat, r, re, El);

    %Plot of Azimuth 
    figure(6+k-1);
    set(gcf,'Menubar','default','Name','Azimuth and Elevation ');
    subplot(1,2,1);
    for j = 1:n_sat
        plot(t, azimuth{j}, color(j));
        hold on;
    end
    plot(t, best_azimuth, 'b', 'LineWidth', 3);
    title('$$Azimuth(t)$$','interpreter','latex');
    xlabel('Time'); % x-axis label
    ylabel('Degree'); % y-axis label
    grid on;

    %Plot of Elevation 
    subplot(1,2,2); 
    for j = 1:n_sat
        plot(t, elevation{j}, color(j));
        hold on;
    end
    plot(t, best_elevation, 'b', 'LineWidth', 3);
    title('$$Elevation(t)$$','interpreter','latex');
    xlabel('Time'); % x-axis label
    ylabel('Degree'); % y-axis label
    grid on;
end

% % % pks = findpeaks(long_deg(1:ceil(length(long_deg)/2))); %pks = most eastern longitude
% % % max_long = min(pks); %the smaller extreme in longitude
% % % locs = find(long_deg == max_long); %index of max_long in the long vector
% % % 
% % % ext_lat = lat_deg(locs); %ext_lat = equivalent latitude at max_long
% % % lat_index = find(lat_deg == ext_lat); %lat_index = index of the most western longitude
% % % if lat_index(1) == locs
% % %     lat_locs = lat_index(2);
% % % else
% % %     lat_locs = lat_index(1);
% % % end
% % % min_long = long_deg(lat_locs); %most western longitude
% % % radius = lambda(locs); %radius of the circle of coverage for the max_long point
% % % rext_east = max_long + radius;
% % % lext_east = min_long + radius;
% % % semi_cov = abs(rext_east - lext_east);
% % % 
% % % cov = 2*radius; %actual double coverage for tundra orbit in the north part
% % % n_orbits = ceil(200/cov);
