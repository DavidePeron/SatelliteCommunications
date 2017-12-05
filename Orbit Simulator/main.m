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
load('molniya.mat');

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
El = 0; %Minimum elevation in degrees
n_sat = 2;
color = ['r','g','c','y','w'];

%Tundra's Raan
% raan = deg2rad(120);
%Molniya's Raan
raan = deg2rad(25);
% lat_gs = [61.267865]; %Ground station in Evenkiysky (Russia)
% long_gs = [-96.608223];

%Create a ground station for each cell (cell of 10 deg of longitude and 10
%of latitude

% long_gs = linspace(-180,0, 19);
% long_gs = repelem(long_gs,3);
% 
% lat_gs = zeros(1,length(long_gs));
% 
% for i = 1:3:length(long_gs)
%     lat_gs(i) = 65;
%     lat_gs(i+1) = 75;
%     lat_gs(i+2) = 85;
%     
% end
% 
% mean_el = zeros(1,length(lat_gs));
lat_gs = [65]; %Ground station in Evenkiysky (Russia)
long_gs = [-90];

R = [22.2016]; %Rain attenuation coefficient [mm/h]
R = repelem(R,length(lat_gs));
lat_us = lat_gs; %user random positioned
long_us = long_gs;
% raan_sat = [raan, raan + deg2rad(-242.4783429000001), raan + deg2rad(-122.1497972999974)];
raan_sat = zeros(1,n_sat)+raan; % Initial raan vector
t = 0:runspeed:T*n_rev; % time [s]

r = cell(n_sat); %cell vector containing the coordinates of the three satellites
phi = cell(n_sat);
best_r = zeros(1,length(t));
best_sat_matrix = zeros(length(lat_gs),length(t));
lambda_matrix = zeros(length(lat_gs),length(t));
lambda_us_matrix = zeros(length(lat_gs),length(t));
CNtot = zeros(length(lat_gs),length(t));
link_margin = zeros(length(lat_gs),length(t));
CNtot_rain = zeros(length(lat_gs),length(t));
link_margin_rain = zeros(length(lat_gs),length(t));

%% RAAN ESTIMATION FOR DIFFERENT SATELLITE IN THE SAME ORBITAL PLANE%

%Computation of polar coordination
[r_t, phi_t] = get_polar(t, T, e, a, mu, n_rev);

r{1} = r_t;
phi{1} = phi_t;

for i=1:n_sat-1    
    r{i+1} = circshift(r_t,(T/(n_sat)*i)./runspeed + 1,2);
    phi{i+1} = circshift(phi_t,(T/(n_sat)*i)./runspeed + 1,2);
end

%In questo modo ogni satellite ha il raan sfasato dello stesso angolo
for i=2:n_sat
    raan_sat(i) = raan_sat(i) + deg2rad(180/n_sat*(i-1));
end

%% COMPUTATION OF THE TRAJECTORY%

[ECI, long, lat] = polar_to_LLA(t, r, phi, inc, w, raan_sat);

%% 3-D ANIMATION

% % %Initializing the Drawing Space and static components
earthmap = imread('planisphere2.jpg');
% 
% earthplot3D(earthmap, t, r, re, El, e, a, ECI, color);
% 
% %% Ground Track
% 
% groundtrack_plot(earthmap, t, r, re, El, long, lat, color, lat_gs, long_gs);

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

%% Azimuth estimation

for k = 1:length(lat_gs)
    
    [azimuth, elevation, best_elevation, best_azimuth, best_sat_matrix(k,:), lambda_matrix(k,:), lambda_us_matrix(k,:)] = view_angles(t, lat_gs(k), long_gs(k), long, lat, r, re, El, lat_us(k),long_us(k));
    
%     mean_el(k) = mean(best_elevation);
    
    [CNtot(k,:),link_margin(k,:), CNtot_rain(k,:),link_margin_rain(k,:)]=link_budget(r,t,lambda_matrix(k,:),lambda_us_matrix(k,:), best_sat_matrix(k,:), R(k), best_elevation);
    %Plot of Azimuth 
    figure(5+k);
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
    
    for j = 1:length(t)
        best_r(j) = r{best_sat_matrix(k,j)}(j);
    end
    
    figure(5+length(lat_gs)+k)
    set(gcf,'Menubar','default','Name','Link Budget');
    subplot(2,2,1);
    plot((best_r./1000),CNtot(k,:))
    title('Distance Best Satellite-Earth VS Total SNR')
    xlabel('Distance Satellite-Earth (Km)')
    ylabel('Total SNR (dB)')
    grid on;

    subplot(2,2,2);
    plot(best_r./1000,link_margin(k,:))
    title('Distance Best Satellite-Earth VS Link margin')
    xlabel('Distance Satellite-Earth (Km)')
    ylabel('Link margin')
    grid on;
    
    subplot(2,2,3);
    plot(t,CNtot(k,:))
    title('Time VS Total SNR')
    xlabel('Time (s)')
    ylabel('Total SNR (dB)')
    grid on;
    
    subplot(2,2,4);
    plot(t,link_margin(k,:))
    title('Time VS Link margin')
    xlabel('Time (s)')
    ylabel('Link margin')
    grid on;
    
    figure(5+2*length(lat_gs)+k)
    set(gcf,'Menubar','default','Name','Link Budget with rain');
    subplot(2,2,1);
    plot((best_r./1000),CNtot_rain(k,:))
    title('Distance Best Satellite-Earth VS Total SNR')
    xlabel('Distance Best Satellite-Earth (Km)')
    ylabel('Total SNR (dB)')
    grid on;

    subplot(2,2,2);
    plot(best_r./1000,link_margin_rain(k,:))
    title('Distance Best Satellite-Earth VS Link margin')
    xlabel('Distance Best Satellite-Earth (Km)')
    ylabel('Link margin')
    grid on;
    
    subplot(2,2,3);
    plot(t,CNtot_rain(k,:))
    title('Distance Time VS Total SNR')
    xlabel('Time (s)')
    ylabel('Total SNR (dB)')
    grid on;
    
    subplot(2,2,4);
    plot(t,link_margin_rain(k,:))
    title('Distance Time VS Link margin')
    xlabel('Time (s)')
    ylabel('Link margin')
    grid on;
    
end

% [value, I] = max(mean_el);
% disp("Best Ground Station at: ");
% latitude_gs = lat_gs(I)
% longitude_gs = long_gs(I)