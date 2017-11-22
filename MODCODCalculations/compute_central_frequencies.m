clear all;
close all;

n_trans = 47;

% Total Downlink Band from 10.7 GHz to 12.7 GHz => 2 GHz available
low_limit_down_band = 10.7e9;
high_limit_down_band = 12.7e9;
down_band = high_limit_down_band - low_limit_down_band;
B_down = 41.6e6;
spatial_division_down = B_down / n_trans - 1;

% Total Uplink Band from 14 GHz to 14.5 GHz => 500 MHz available
low_limit_up_band = 14e9;
high_limit_up_band = 14.5e9;
up_band = high_limit_up_band - low_limit_up_band; 
B_up = 10.4e6;
spatial_division_up = B_up / n_trans - 1;

centrals_down = zeros(1,n_trans - 1);
for i = 1:n_trans
    centrals_down(i) = B_down/2 + (B_down + B_down/(n_trans - 1))*(i-1);
end

centrals_down = centrals_down + low_limit_down_band;

centrals_up = zeros(1,n_trans - 1);
for i = 1:n_trans
    centrals_up(i) = B_up/2 + (B_up + B_up/(n_trans - 1))*(i-1);
end

centrals_up = centrals_up + low_limit_up_band;