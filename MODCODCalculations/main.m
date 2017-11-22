clear all;
close all;

% 52 Mhz Transponders
B_trans = 52e6;

% Total Downlink Band from 10.7 GHz to 12.7 GHz => 2 GHz available
low_limit_down_band = 10.7e9;
high_limit_down_band = 12.7e9;
down_band = high_limit_down_band - low_limit_down_band;

% Total Uplink Band from 14 GHz to 14.5 GHz => 500 MHz available
low_limit_up_band = 14e9;
high_limit_up_band = 14.5e9;
up_band = high_limit_up_band - low_limit_up_band;

rolloff_factor = 0.15;

%Bit rate Downlink and Uplink in bit per second
bit_rate_down = 50e6;
bit_rate_up = 5e6;

max_efficiency = 5.900855;

% % 36 Mhz Transponders
% n_trans = 70;
% B_trans = 36e6;

nu = 0.5;
step = 0.02;
epsilon = 2e5;

contained = false;
above = false;
converge = false;

%% Find the first solution, with a lot of black space in the available
% bandwidth
while(~contained)
    Rsymb_downlink = bit_rate_down/nu;
    B_downlink = Rsymb_downlink*(1+rolloff_factor);
    if(B_downlink < B_trans)
        Rsymb_uplink = bit_rate_up/nu;
        B_uplink = Rsymb_uplink*(1+rolloff_factor);
        if(B_uplink + B_downlink < B_trans)
            contained = true;
        end
    end
    nu = nu + step;
end

space = B_trans - B_uplink - B_downlink;
space_old = space;
%% Try to fulfill the blackspace of available bandwidth decreasing the spectral efficiency (nu)
step = step/2;
while(~converge)
    while(~above)
        nu = nu - step;
        Rsymb_downlink = bit_rate_down/nu;
        B_downlink = Rsymb_downlink*(1+rolloff_factor);
        if(B_downlink < B_trans)
            Rsymb_uplink = bit_rate_up/nu;
            B_uplink = Rsymb_uplink*(1+rolloff_factor);
            space_old = space;
            space = B_trans - B_uplink - B_downlink;
            if(space < 0)
                above = true;
                space = space_old;
            end
        end
    end
    nu = nu + step;
    if(space < epsilon)
        converge = true;
    end
end
 
Rsymb_downlink = bit_rate_down/nu;
B_downlink = Rsymb_downlink*(1+rolloff_factor);
Rsymb_uplink = bit_rate_up/nu;
B_uplink = Rsymb_uplink*(1+rolloff_factor);

%% Calculate the number of transponders needed
n_trans = floor(down_band/B_downlink);

centrals = zeros(1,n_trans - 1);
for i = 1:n_trans - 1
    centrals(i) = B_downlink / 2 + (B_downlink + B_downlink/(n_trans - 2))*(i-1);
end

centrals = centrals + low_limit_down_band;

%% Calculation of the actual spectral efficiency 

max_carriers = 8;
n_carriers = linspace(1, max_carriers, max_carriers);
efficiency = zeros(1,max_carriers);
Band_for_carrier_down = zeros(1,max_carriers);
Band_for_carrier_up = zeros(1,max_carriers);

for i = 1:max_carriers
    Band_for_carrier_down(i) = B_downlink/i;
    Band_for_carrier_up(i) = B_uplink/i;
    efficiency(i) = (1+rolloff_factor)*bit_rate_down/Band_for_carrier_down(i);
end

figure(1);
stem(n_carriers, efficiency);
grid on;
title("Efficiency considering different number of carriers");
xlabel("$$\# ~Carriers$$", 'Interpreter' , 'latex');
ylabel("Spectral efficiency $$\nu (bit/s/Hz)$$", 'Interpreter' , 'latex');

%Find the best suitable number of carriers
suitable_num_carriers_index = efficiency < max_efficiency;
suitable_num_carriers = efficiency(suitable_num_carriers_index);
best_efficiency = suitable_num_carriers(end);
best_n_carriers = length(suitable_num_carriers);

%Now we use the best efficiency and we take the modulation in the
%DVB-S2X table that has the nearest value of spectral efficiency to the
%best one (the greatest one)