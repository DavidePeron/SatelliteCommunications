clear all;
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
El = 0; %Minimum elevation in degrees
n_sat = 2;
color = ['r','g','c','y','w'];
raan = deg2rad(120);
% raan_sat = [raan, raan + deg2rad(-242.4783429000001), raan + deg2rad(-122.1497972999974)];
raan_sat = zeros(1,n_sat)+raan; % Initial raan vector
t = 0:runspeed:T*n_rev; % time [s]

r = cell(n_sat); %cell vector containing the coordinates of the three satellites
phi = cell(n_sat);
h1 = cell(n_sat);

%% RAAN ESTIMATION FOR DIFFERENT SATELLITE IN THE SAME ORBITAL PLANE%

%Computation of polar coordination
[r_t, phi_t] = get_polar(t, T, e, a, mu, n_rev);

for i=0:n_sat-1    
    r{i+1} = circshift(r_t,(T/n_sat*i)./runspeed,2);
    phi{i+1} = circshift(phi_t,(T/n_sat*i)./runspeed,2);
    h1{i+1} = plot(0,0,'o','Markersize',1);
end

[~, long, ~] = polar_to_LLA(t, r, phi, inc, w, raan_sat);

for i=2:n_sat
    long_temp = circshift(long{1},(T/n_sat*(i-1))./runspeed,2);
    raan_diff = long{1}(1) - long{i}(1);
    raan_sat(i) = raan_sat(i) + deg2rad(raan_diff);
end

%% COMPUTATION OF THE TRAJECTORY%

[ECI, long, lat] = polar_to_LLA(t, r, phi, inc, w, raan_sat);

%% 3-D ANIMATION

%Initializing the Drawing Space and static components
earthmap = imread('planisphere2.jpg');
set(gcf,'Menubar','default','Name','Orbit Visualization', ... 
    'NumberTitle','off','Position',[70,10,750,750]); 
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
[xx yy zz]=ellipsoid (0,0,0,a_WGS84, a_WGS84, b);
earthmap2 = flip(earthmap,1);
earthmap2 = flip(earthmap2,2);
pro.FaceColor= 'texture';
pro.EdgeColor = 'none';
pro.FaceLighting = 'phong';
pro.Cdata = earthmap2;
earth= surface(xx,yy,zz,pro);

rotate (earth, [0 0 1], 0);
Xaxis= line([-1e7 1e7],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis= line([0 0],[-1e7 1e7],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
rotate (Xaxis, [0 0 1], 0);
rotate (Yaxis, [0 0 1], 0);
% Sun=light('Position',[1 0 0],'Style','infinite');

%Plotting the ECI Axes
line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')

sphere_position = cell(n_sat);
position = cell(n_sat);

for i=1:length(t)
    rotate (earth, [0 0 1], (360/length(t)), [0,0,0]);
    rotate (Xaxis, [0 0 1], (360/length(t)), [0,0,0]);
    rotate (Yaxis, [0 0 1], (360/length(t)), [0,0,0]);

    %Drawing the red sphere
    for j=1:n_sat
        %Drawing the red sphere
        sphere_position{j}(i)=plot3 (ECI{j}(1,i), ECI{j}(2,i), ECI{j}(3,i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',color(j),'MarkerSize', 6);
        position{j}(i)=line([0 ECI{j}(1,i)],[0 ECI{j}(2,i)], [0 ECI{j}(3,i)],'Color', color(j), 'LineWidth', 2);
        
        if(i~=1)
            set (sphere_position{j}(i-1), 'Visible', 'off');
            set (position{j}(i-1), 'Visible', 'off');
        end
        
        if (i~=1 && i<=length(t))
            line([ECI{j}(1,i-1) ECI{j}(1,i)],[ECI{j}(2,i-1) ECI{j}(2,i)], [ECI{j}(3,i-1) ECI{j}(3,i)],'Color', 'black', 'LineWidth', 1.5);
        end
    end
    
    pause (0.01);
end

%% Ground Track
figure (2);
set(gcf,'Menubar','none','Name','Earth Track', ... 
    'NumberTitle','off','Position',[70,30,1000,500]); 
hold on
image([180 -180],[90 -90],earthmap,'CDataMapping', 'scaled');
axis equal
axis ([-180 180 -90 90]);
set(gca, 'XDir', 'reverse');

%% Ground Stations plot

% Earth-satellite geometry parameters
elong_ext = zeros(1,length(t)); %longitude of the extreme toward east
nlat_ext = zeros(1,length(t)); %latitude of the extreme toward north
wlong_ext = zeros(1,length(t)); %longitude of the extreme toward west
slat_ext = zeros(1,length(t)); %latitude of the extreme toward south

% plot (167.717,8.717,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
% plot (360-76.496, 42.440,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
for i=1:length(t)
    for j = 1:n_sat
        %% Nadir angle estimation
        rho = asin(re/r{j}(i));
        beta = asin(cos(deg2rad(El))*sin(rho));
        nadir = beta;

        %% Lambda and arch estimation 
        lambda= 90 - rad2deg(nadir) - El; %central angle in degrees

        arch = (2*pi*r{j}(i)/360)*lambda; %arch between subsat. point and extreme position of a user
        l = 2*re*sin(deg2rad(lambda/2)); %segment between subsat. point and extreme position of a user

        %Finding the LLA coordinates of the extreme position of a user
        %We consider positive angles towards east and north and negative ones 
        %towards west and south; 
        elong_ext = (long{j}(i) + lambda) - 360*floor((long{j}(i) + lambda)/360);
        wlong_ext = (long{j}(i) - lambda);
        nlat_ext = (lat{j}(i) + lambda); %- 90*floor(lat{j}(i)+lambda/90);
        slat_ext = lat{j}(i) - lambda;
    
        plot (long{j}(i),lat{j}(i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',color(j),'MarkerSize', 2);
        if (i~=1 && abs(long{j}(i-1)-long{j}(i))<100)
            line([long{j}(i-1) long{j}(i)],[lat{j}(i-1) lat{j}(i)],'Color', color(j), 'LineWidth', 2);
        end
        delete(h1{j});
        h1{j} = plot(long{j}(i),lat{j}(i),'o','MarkerEdgeColor',color(j), 'Markersize', lambda*2);
    end
    
    pause (0.01);
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
% 
% %Plot of Anomalies
% figure(7);
% set(gcf,'Menubar','default','Name','Anomalies'); 
% plot(t, M, 'r');
% hold on
% plot(t, E, 'b');
% hold on
% plot(t, phi_t, 'g');

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
prompt = '\n Insert a latitude for the ground station: ';
lat_gs = input(prompt);
prompt = '\n Insert a longitude between [-180; +180] for the ground station: ';
long_gs = input(prompt);
azimuth = cell(n_sat);
a_primo = zeros(1,length(t));
elevation = cell(n_sat);

for i = 1:length(t)
    for j = 1:n_sat
        cos_lambda = cosd(lat_gs)*cosd(lat{j}(i))*cosd(long_gs - long{j}(i)) + sind(lat_gs)*sind(lat{j}(i));
        lambda = acosd(cos_lambda);
        a_az = asind(sind(abs(long_gs - long{j}(i)))*cosd(lat_gs)/sind(lambda));

        if lat{j}(i) < lat_gs && long{j}(i) < long_gs % South east
            azimuth{j}(i) = 180 - a_az;
        end
        if lat{j}(i) < lat_gs && long{j}(i) > long_gs % South west
            azimuth{j}(i) = 180 + a_az;
        end
        if lat{j}(i) > lat_gs && long{j}(i) < long_gs % North east
            azimuth{j}(i) = a_az;
        end
        if lat{j}(i) > lat_gs && long{j}(i) > long_gs % North west
            azimuth{j}(i) = 360 - a_az;
        end    
%         a_primo(i) = a_az;

        b = 1+(re/r{j}(i))^2;
        c = 2*re/r{j}(i);
        cos_gamma = (c*(cosd(El))^2+sqrt(c^2*(cosd(El))^4-4*b*(cosd(El))^2+4))/2;
        gamma = acosd(cos_gamma);
        if(lambda < gamma)
            elevation{j}(i) = rad2deg(acos(sin(acos(cos_lambda))/sqrt(1+(re/r{j}(i))^2+2*re/r{j}(i)*cos_lambda)));
        else
            elevation{j}(i) = 0;
        end
    end
end

max_el = zeros(1,length(t));
temp = zeros(1,n_sat);

for i = 1:length(t)
    for j=1:n_sat
        temp(j) = elevation{j}(i);
    end
    max_el(i) = max(temp);
end

%Plot of Azimuth 
figure(6);
set(gcf,'Menubar','default','Name','Azimuth and Elevation');
subplot(1,2,1);
for j = 1:n_sat
    plot(t, azimuth{j}, color(j));
    hold on;
end
hold off;
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
plot(t, max_el, 'b', 'LineWidth', 3);
title('$$Elevation(t)$$','interpreter','latex');
xlabel('Time'); % x-axis label
ylabel('Degree'); % y-axis label
grid on;

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
