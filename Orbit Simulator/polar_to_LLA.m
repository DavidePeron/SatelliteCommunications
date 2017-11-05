function [ECI, long, lat] = polar_to_LLA(t, r, phi, inc, w, raan_sat)

    load('earth_constants.mat');
    n_sat = length(raan_sat);
    ECI = cell(n_sat);
    long = cell(n_sat);
    lat = cell(n_sat);
    ECEF_coord = zeros(3,length(t));
    x_0 = zeros(1,length(t));
    y_0 = zeros(1,length(t));
    z_0 = zeros(1,length(t));

    %Julian Date of 10 January 2017 at noon
    A = 2017; %Year
    DTA = 10; %Days from 1 Jan of year A
    NAB1900 = 29; %Number of leap years since 1900
    TU = 12; %Hours in the day
    JD = 2415020 + 365*(A - 1900) + DTA + NAB1900 + TU/24 - 0.5;

    T_c = (JD - 2415020)/36525;
    alpha_go = 99.6909833 + 36000.7689*T_c + 3.8708e-4*T_c^2;
    
    %ECEF coordinates and LLA
    for i = 1:length(t)
        %Transformation Matrix from ECI to ECEF Coordinates
        theta = alpha_go + 0.25068447*t(i)/60;
        theta = wrapTo2Pi(deg2rad(theta));
        ECI_to_ECEF = [cos(theta) sin(theta) 0; 
                      -sin(theta) cos(theta) 0; 
                           0          0      1];
        %Compute ECEF Coordinates
        for j = 1:n_sat
            x_0(i) = r{j}(i)* cos(phi{j}(i));
            y_0(i) = r{j}(i)* sin(phi{j}(i));
            
            %Transformation Matrix from orbital to ECI Coordinates
            orbital_to_ECI = [cos(raan_sat(j))*cos(w) - sin(raan_sat(j))*cos(inc)*sin(w), -cos(raan_sat(j))*sin(w) - sin(raan_sat(j))*cos(inc)*cos(w), sin(raan_sat(j))*sin(inc);
            sin(raan_sat(j))*cos(w) + cos(raan_sat(j))*cos(inc)*sin(w), -sin(raan_sat(j))*sin(w) + cos(raan_sat(j))*cos(inc)*cos(w), -cos(raan_sat(j))*sin(inc);
            sin(inc)*sin(w), sin(inc)*cos(w), cos(inc)];


            ECI{j}(:,i) = orbital_to_ECI*[x_0(i); y_0(i); z_0(i)];

            ECEF_coord(:,i) = ECI_to_ECEF*ECI{j}(:,i);

            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transformation from ECEF to LLA
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p = sqrt(ECEF_coord(1,i).^2 + ECEF_coord(2,i).^2);
            th = atan2(a_WGS84*ECEF_coord(3,i), b*p); %theta
            long{j}(i) = atan2(ECEF_coord(2,i), ECEF_coord(1,i));
            lat{j}(i) = atan2(ECEF_coord(3,i)+ep^2.*b.*sin(th).^3 , p-e_WGS84^2.*a_WGS84.*cos(th).^3);

            N = a_WGS84./sqrt(1-e_WGS84.^2*sin(lat{j}(i)).^2);

            long{j}(i) = rad2deg(long{j}(i));
            lat{j}(i) = rad2deg(lat{j}(i));

            %return long in [0,360] range and remove 360 to longitude higher then
            %180, to have West coordinates
             long{j}(i) = mod(long{j}(i),360);
             if long{j}(i) > 180
                 long{j}(i) = long{j}(i) - 360;
             end
        end
    end

end