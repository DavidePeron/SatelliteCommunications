%In this case, we consider the link between GT to the user. It means the
%forward link takes into account the uplink and downliink part.
    
function [CNtot, link_margin, CNtot_rain, link_margin_rain]= link_budget(r,t,lambda,lambda_us, best_sat, R, best_elevation)

    %Vectors
    distance = zeros(1,length(t)); %Distance [Km] between the earth in each instant of time and the Earth (GroundStation position)
    distance_down = zeros(1,length(t));
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
    diameter_sat = 1.2; %Satellite Antenna Diameter [m]
    Diameter_user = 1; %User Antenna Diameter [m]
    diameter_gs = 6; %Ground station Antenna Diameter [m]
    T_user_antenna = 80; %User antenna temperature [K]
    T_antenna_sat = 290; %Sat antenna temperature
    pointing_loss = 0.3; %[dB]
    gases_absortion = 0.3;
    re = 6371; %ground station to earth centre (Km)
    c = 3e8;
    IBO = -16; %[dB]
    OBO = IBO + 6 - 6*exp(IBO/6);
    l_mc = -OBO; %Losses due to multicarrier operation [dB]
    l_ftx = 0.5; %Losses due to feeder [dB]
    carrier_b_r = 50e6; %[bit/s]
    G_lna = 30; %[dB]
    g_lna = 10^(G_lna/10);
    f_lna = 4; %[dB]
    t_lna = 290*((10^(f_lna/10))-1);
    f_cable = 3; %[dB]
    g_cable = 1/10^(f_cable/10);
    t_cable = 290*((10^(f_cable/10))-1);
    f_dc = 10; %[dB] Noise Figure of the down converter
    t_dc = 290*((10^(f_dc/10))-1);
    
    %% Interferences
    carr_to_inter_noise = 23;
    add_deg = 0.5;
    imp_margin = 1.2;
    required_EsNo=6.55;
    %losses = 0.2; %[dB]
    %SFD = -92.5; %[dBW/m^2]
    
    %% Distance calculation GW to satellite and satellite to user
    for j = 1:length(t)   
       distance(j) = r{best_sat(j)}(j)*sqrt(1 + (re/((r{best_sat(j)}(j)./1000))).^2 - 2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda(j)));
       distance_down(j) = ((r{best_sat(j)}(j))*sqrt(1 + (re/((r{best_sat(j)}(j)./1000))).^2 - 2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda_us(j))));
    end
    
    %% UPLINK    
    mod_eff = 2.10485;
    Transp_BW = 72; %size of one transponder [MHz]
    k = 10*log10(1.38e-23); %Boltzmann constant
    frequency_uplink = 14414e6; %Central frequency for uplink in Ku band [MHz]
    p_HPA = 100; %Power of a single HPA multipied by the n of channels [W]
    efficiency = 0.6; %Antenna efficiency 0.6-0.7
    
    carrier_s_r = carrier_b_r/4; %Carrier symbol rate: with 16APSK we have 4 bits/symbol 
    lambda_link_up = c/(frequency_uplink); %Uplink wavelength 
    Ggw = 10*log10(efficiency*(pi*diameter_gs/lambda_link_up)^2); %Gain of the Ground station antenna [dBi]
    p_t = 20*log10(p_HPA)-l_mc-l_ftx; %Total power for a transponder
    p_tx = p_t + 10*log10(12); %Total power tx by the antenna (p_t x the n of transponders)
    G_rx_sat = 10*log10(efficiency*(pi*diameter_sat/lambda_link_up)^2); %Gain of the Satellite antenna in reception [dBi]
    T_system_sat = T_antenna_sat + t_lna + t_cable/g_lna + t_dc/(g_cable*g_lna);
    satellite_G_T = G_rx_sat - 10*log10(T_system_sat); %G/T satellite
    EIRPgw = Ggw + p_tx; %Gateway EIRP
    
    for i=1:length(t) 
        Lu_pathloss(i) = 20*log10(4*pi*distance(i)*frequency_uplink/c);
        Lu_rain(i) = (4.21e-5*(frequency_uplink/10^9)^2.42)*R^(1.41*(frequency_uplink/10^9)^(-0.0779))*2/sind(best_elevation(i)); %I assume that the clouds are 2km high
        CN_uplink(i) = EIRPgw - Lu_pathloss(i) - k - 10*log10(carrier_s_r) + satellite_G_T;
        CN_uplink_rain(i) = CN_uplink(i) - Lu_rain(i);
    end

    %% DOWNLINK
    frequency_downlink = 11114e6;%Higher central frequency for downlink in Ku band
    p_HPA_down = 20; %Power of a single HPA multiplied by the n of channels
    
    lambda_link_dw = c/frequency_downlink; %Downlink wavelength
    G__tx_sat = 10*log10(efficiency*(pi*diameter_sat/lambda_link_dw)^2); %Gain of the Satellite antenna in tx [dBi]
    p_t_down = 20*log10(p_HPA_down) - l_mc - l_ftx; %Total power for a transponder
    p_tx_down = p_t_down + 10*log10(12); %Total power tx by the antenna (p_t x the n of transponders)
    EIRP_down = G__tx_sat + p_tx_down; %Gateway EIRP 
    G_user = 10*log10(efficiency*((pi*Diameter_user*frequency_downlink)/c)^2);
    T_system_gs = T_user_antenna + t_lna + t_cable/g_lna + t_dc/(g_cable*g_lna);
    GT_user = G_user - 10*log10(T_system_gs); %G/T user station

    for n = 1:length(t)   
        Ld_pathloss(n) = 20*log10((4*pi*distance_down(n)*frequency_downlink)/c);  %Path Loss
        Ld_rain(n) = (4.21e-5*(frequency_downlink/10^9)^2.42)*R^(1.41*(frequency_downlink/10^9)^(-0.0779))*2/sind(best_elevation(n)); %I assume that the clouds are 2km high
        CN_downlink(n) = EIRP_down + GT_user - Ld_pathloss(n) - k - 10*log10(carrier_s_r) - pointing_loss - gases_absortion;
        CN_downlink_rain(n) = CN_downlink(n) - Ld_rain(n);
    end

    %% CN_TOT    
    for m = 1:length(t)
        CNtot(m) = -20*log10(1/(10^(CN_uplink(m)/20)) + 1/(10^(CN_downlink(m)/20) + 1/(10^(carr_to_inter_noise/20))));
        available_EsNo(m) = CNtot(m) - add_deg - imp_margin;
        link_margin(m) = available_EsNo(m) - required_EsNo;
        CNtot_rain(m) = -20*log10(1/(10^(CN_uplink_rain(m)/20)) + 1/(10^(CN_downlink_rain(m)/20) + 1/(10^(carr_to_inter_noise/20))));
        available_EsNo_rain(m) = CNtot_rain(m) - add_deg - imp_margin;
        link_margin_rain(m) = available_EsNo_rain(m) - required_EsNo;
    end
end