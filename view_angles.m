function [azimuth, elevation, best_elevation, best_azimuth] = view_angles(t, lat_gs, long_gs, long, lat, r, re, El)

    n_sat = length(long);
    azimuth = cell(n_sat);
    elevation = cell(n_sat);
    temp = zeros(1,n_sat);
    best_elevation = zeros(1,length(t));
    best_azimuth = zeros(1,length(t));

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

            b = 1+(re/r{j}(i))^2;
            c = 2*re/r{j}(i);
            cos_gamma = (c*(cosd(El))^2+sqrt(c^2*(cosd(El))^4-4*b*(cosd(El))^2+4))/2;
            gamma = acosd(cos_gamma);
%             elevation{j}(i) = rad2deg(acos(sin(acos(cos_lambda))/sqrt(1+(re/r{j}(i))^2+2*re/r{j}(i)*cos_lambda)));
            if(lambda < gamma)
                elevation{j}(i) = rad2deg(acos(sin(acos(cos_lambda))/sqrt(1+(re/r{j}(i))^2+2*re/r{j}(i)*cos_lambda)));
            else
                elevation{j}(i) = 0;
            end
            
            temp(j) = elevation{j}(i);
        end
        
        % take the higher elevation between the different satellites
        [best_elevation(i), best_sat] = max(temp);
        best_azimuth(i) = azimuth{best_sat}(i);
    end

end