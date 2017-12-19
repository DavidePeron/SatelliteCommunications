%In this case, we consider the link between GT to the user. It means the
%forward link takes into account the uplink and downliink part.
%%%%%
%INPUT PARAMETERS:
%r = distance between each satellite and the Earth center;
%t = time vector;
%lambda = central angle of the gs;
%lambda_us = central angle of the user station;
%best_sat = vector containing which of the two satellite has the best
%elevation for each instant of time;
%R = rain factor for 0.01 prob.;
%best_elevation = elevation value corresponding to the best_sat element;
%%%%%    

function [CNtot, link_margin, CNtot_rain, link_margin_rain]= link_budget(r,t,lambda,lambda_us, best_sat, R_up, R_down, best_elevation)

    %Vectors
    distance = zeros(1,length(t)); %Distance [km] between the satellite in each instant of time and the Earth (GroundStation position)
    distance_down = zeros(1,length(t)); %Distance of the user station [km]
    Lu_pathloss = zeros(1,length(t)); %Uplink Losses
    Lu_rain = zeros(1,length(t)); %Uplink Rain Losses
    Ld_pathloss = zeros(1,length(t)); %Downlink Losses
    Ld_rain = zeros(1,length(t)); %Downlink Rain Losses
    CN_uplink = zeros(1,length(t)); %SNR Uplink
    CN_downlink = zeros(1,length(t)); %SNR Downlink
    CN_uplink_rain = zeros(1,length(t)); %SNR Uplink with rain
    CN_downlink_rain = zeros(1,length(t)); %SNR Downlink with rain
    CNtot = zeros(1,length(t)); %Total SNR of the link budget
    available_EsNo = zeros(1,length(t));
    link_margin = zeros(1,length(t));
    CNtot_rain = zeros(1,length(t)); %Total SNR of the link budget
    available_EsNo_rain = zeros(1,length(t));
    link_margin_rain = zeros(1,length(t));
    
    %% Constants
    diameter_sat = 1; %Satellite Antenna Diameter [m]
    Diameter_user = 1; %User Antenna Diameter [m]
    diameter_gs = 3; %Ground station Antenna Diameter [m]
    T_user_antenna = 80; %User antenna temperature [K]
    T_antenna_sat = 290; %Sat antenna temperature
    pointing_loss = 0.3; %[dB]
    gases_absortion = 0.3;
    re = 6371; %ground station to earth centre (Km)
    c = 3e8;
    IBO = -0.5; %[dB]
    OBO = IBO + 6 - 6*exp(IBO/6);
    l_mc = -OBO; %Losses due to multicarrier operation [dB]
    l_ftx = 0.5; %Losses due to feeder [dB]
    carrier_b_r = 50e6; %carrier bit rate [bit/s]
    G_lna = 30; %LNA gain [dB]
    g_lna = 10^(G_lna/10); %LNA gain
    f_lna = 4; %LNA noise figure [dB]
    t_lna = 290*((10^(f_lna/10))-1); %LNA noise temperature
    f_dc = 10; %[dB] Noise Figure of the down converter
    t_dc = 290*((10^(f_dc/10))-1);%down converter temperature
    
    %% Interferences
    carr_to_inter_noise = 23; %intermodulation noise
    add_deg = 0.5; %addicitonal degradation
    imp_margin = 1.2; %implementation margin
    required_EsNo = 6.55; %Taken from table 20a of DVB-S2x standard document 
    %(We are not sure about the correctness of the selected value since it does not consider the role
    % of the interference but it assumes an AWGN channel)
    
    %% Distance calculation GW to satellite and satellite to user
    for j = 1:length(t)   
       distance(j) = r{best_sat(j)}(j)*sqrt(1 + (re/((r{best_sat(j)}(j)./1000))).^2 - 2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda(j)));
       distance_down(j) = ((r{best_sat(j)}(j))*sqrt(1 + (re/((r{best_sat(j)}(j)./1000))).^2 - 2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda_us(j))));
    end
    
    %% UPLINK    
    k = 10*log10(1.38e-23); %Boltzmann constant
    frequency_uplink = 14414e6; %Central frequency for uplink in Ku band [MHz]
    p_HPA = 100; %Power of a single HPA multipied by the n of channels [W]
    efficiency = 0.6; %Antenna efficiency 0.6-0.7
    
    carrier_s_r = carrier_b_r/4; %Carrier symbol rate: with 16APSK we have 4 bits/symbol 
    lambda_link_up = c/(frequency_uplink); %Uplink wavelength 
    Ggw = 10*log10(efficiency*(pi*diameter_gs/lambda_link_up)^2); %Gain of the Ground station antenna [dBi]
    p_tx = 10*log10(p_HPA)-l_mc-l_ftx; %Total power tx by antenna
    G_rx_sat = 10*log10(efficiency*(pi*diameter_sat/lambda_link_up)^2); %Gain of the Satellite antenna in reception [dBi]
    T_system_sat = T_antenna_sat + t_lna + t_dc/g_lna;
    satellite_G_T = G_rx_sat - 10*log10(T_system_sat); %G/T satellite
    EIRPgw = Ggw + p_tx; %Gateway EIRP
    
    for i=1:length(t) 
        Lu_pathloss(i) = 20*log10(4*pi*distance(i)*frequency_uplink/c);
        L_s = 2/sind(best_elevation(i)); %I assume that the clouds are 2km high
        r_p = 90/(90+4*L_s*cosd(best_elevation(i))); % Non homogeneity of the rain
        Lu_rain(i) = (4.21e-5*(frequency_uplink/10^9)^2.42)*R_up^(1.41*(frequency_uplink/10^9)^(-0.0779))*L_s*r_p; %I assume that the clouds are 2km high
%         Lu_rain(i) = Lu_rain(i)*(1/0.01)^(-0.5);
        CN_uplink(i) = EIRPgw - Lu_pathloss(i) - k - 10*log10(carrier_s_r) - pointing_loss - gases_absortion + satellite_G_T - abs(IBO);
        CN_uplink_rain(i) = CN_uplink(i) - Lu_rain(i);
    end

    %% DOWNLINK
    frequency_downlink = 11114e6;%Higher central frequency for downlink in Ku band
    p_HPA_down = 20; %Power of a single HPA multiplied by the n of channels
    
    lambda_link_dw = c/frequency_downlink; %Downlink wavelength
    G__tx_sat = 10*log10(efficiency*(pi*diameter_sat/lambda_link_dw)^2); %Gain of the Satellite antenna in tx [dBi]
    p_tx_down = 10*log10(p_HPA_down) - l_mc - l_ftx; %Total power tx by the antenna (p_t x the n of transponders)
    EIRP_down = G__tx_sat + p_tx_down; %Gateway EIRP 
    G_user = 10*log10(efficiency*((pi*Diameter_user*frequency_downlink)/c)^2);
    T_system_gs = T_user_antenna + t_lna + t_dc/g_lna;
    GT_user = G_user - 10*log10(T_system_gs); %G/T user station

    for n = 1:length(t)   
        Ld_pathloss(n) = 20*log10((4*pi*distance_down(n)*frequency_downlink)/c);  %Path Loss
        L_s = 2/sind(best_elevation(n)); %I assume that the clouds are 2km high
        r_p = 90/(90+4*L_s*cosd(best_elevation(n))); % Non homogeneity of the rain
        Ld_rain(n) = (4.21e-5*(frequency_downlink/10^9)^2.42)*R_down^(1.41*(frequency_downlink/10^9)^(-0.0779))*L_s*r_p; %I assume that the clouds are 2km high
        CN_downlink(n) = EIRP_down + GT_user - Ld_pathloss(n) - k - 10*log10(carrier_s_r) - pointing_loss - gases_absortion - abs(OBO);
        CN_downlink_rain(n) = CN_downlink(n) - Ld_rain(n);
    end

    %% CN_TOT    
    for m = 1:length(t)
        CNtot(m) = -10*log10(1/(10^(CN_uplink(m)/10)) + 1/(10^(CN_downlink(m)/10) + 1/(10^(carr_to_inter_noise/10))));
        available_EsNo(m) = CNtot(m) - add_deg - imp_margin;
        link_margin(m) = available_EsNo(m) - required_EsNo;
        CNtot_rain(m) = -20*log10(1/(10^(CN_uplink_rain(m)/20)) + 1/(10^(CN_downlink_rain(m)/20) + 1/(10^(carr_to_inter_noise/20))));
        available_EsNo_rain(m) = CNtot_rain(m) - add_deg - imp_margin;
        link_margin_rain(m) = available_EsNo_rain(m) - required_EsNo;
    end
end