function [CNtot, link_margin, CNtot_rain, link_margin_rain]= link_budget(r,t,lambda,lambda_us, best_sat, R, best_elevation)

    %In this case, we consider the link between GT to the user. It means the
    %forward link taking into account the uplink and downliink part.

    %vectors
    distance = zeros(1,length(t));%Distance [Km] between the earth in each instant of time and the Earth (GroundStation position)
    distance_down = zeros(1,length(t));
    Lu_pathloss = zeros(1,length(t));%Uplink Losses
    Lu_rain = zeros(1,length(t));%Uplink Rain Losses
    Ld_pathloss = zeros(1,length(t));%Downlink Losses
    Ld_rain = zeros(1,length(t));%Downlink Rain Losses
    CN_uplink = zeros(1,length(t));%SNR Uplink
    CN_downlink = zeros(1,length(t));%SNR Downlink
    CN_uplink_rain = zeros(1,length(t));%SNR Uplink with rain
    CN_downlink_rain = zeros(1,length(t));%SNR Downlink with rain
    CNtot = zeros(1,length(t));%Total SNR of the link budget
    available_EsNo=zeros(1,length(t));%
    link_margin=zeros(1,length(t));%
    CNtot_rain = zeros(1,length(t));%Total SNR of the link budget
    available_EsNo_rain=zeros(1,length(t));%
    link_margin_rain=zeros(1,length(t));%
    
    % Distance calculation GW to satellite and satellite to user
    re=6371;%ground station to earth centre (Km)
    for j=1:length(t)   
       distance(j)=r{best_sat(j)}(j)*sqrt(1+(re/((r{best_sat(j)}(j)./1000))).^2-2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda(j)));
       distance_down(j)=((r{best_sat(j)}(j))*sqrt(1+(re/((r{best_sat(j)}(j)./1000))).^2-2*(re/((r{best_sat(j)}(j)./1000)))*cosd(lambda_us(j))));
    end

    %Modulation 16APSK 8/15-L
    %Efficiency 2.10485
    %Es/No required 6.55 dB
    % Transponder 72 MHz
    %UPLINK
    mod_eff = 2.10485;
    carrier_b_r = 50e6; %[bit/s]
    carrier_s_r = carrier_b_r/4; %Carrier symbol rate: with 16APSK we have 4 bits/symbol 
    % carrier_BW = ((carrier_b_r)/(mod_eff))*(1+0.15);% Carrier BW
    Transp_BW = 72; %size of one transponder [MHz]
    k = 10*log10(1.38e-23);%Boltzmann constant
    %BW = 10*log10(carrier_BW); %Carrier Bw in dB
    %losses = 0.2; %[dB]



    % SFD = -92.5; %[dBW/m^2]
    c = 3e8;
    frequency_uplink=14414e6;%Central frequency for uplink in Ku band [MHz]
    lambda_link_up = c/(frequency_uplink);%Uplink wavelength 
    Diameter = 11; %Antenna Diameter
    min_beam_width =70*lambda_link_up/Diameter;
    efficiency = 0.6;%Antenna efficiency 0.6-0.7
    Ggw=efficiency*48360/min_beam_width^2;
    Ggw=10*log10(Ggw);
    IBO=-0.5;
    OBO=IBO+6-6*exp(IBO/6);
    l_mc = -OBO; %Losses due to multicarrier operation [dB]
    l_ftx = 0.5; %Losses due to feeder [dB]
    p_HPA = 100; %Power of a single HPA multipied by the n of channels
    p_t = 20*log10(p_HPA)-l_mc-l_ftx; %Total power for a transponder
    p_tx = p_t + 10*log10(12); %Total power tx by the antenna (p_t x the n of transponders)
    satellite_G_T=6; %[dB]
    EIRPgw = Ggw + p_tx; %Gateway EIRP
    EIRP_mantain_o_p = zeros(1,length(t));
    Required_amp_pow = zeros(1,length(t));

    for i=1:length(t)
    %    EIRP_mantain_o_p(i)= SFD + 10*log10(4*pi*(distance(i)).^2) + IBO; 
       Required_amp_pow(i)=10^((EIRP_mantain_o_p(i)-Ggw)/10);
       Lu_pathloss(i)= 20*log10((4*pi*distance(i)*frequency_uplink)/c);
       Lu_rain(i) = (4.21e-5 * (frequency_uplink/10^9)^2.42)*R^(1.41*(frequency_uplink/10^9)^(-0.0779))*2/sind(best_elevation(i)); %I assume that the clouds are 2km high
       %92.44 + 20*log10(frequency_uplink)+20*log10(distance(1,m))+losses; %freq in GHz, d in km, the distance depends on the time. It would be r_t r_t(1,m)./1000
       CN_uplink(i) = EIRPgw -Lu_pathloss(i)-k-10*log10(carrier_s_r)+satellite_G_T;
       CN_uplink_rain(i) = CN_uplink(i) - Lu_rain(i);
    end

    %DOWNLINK

    % carrier_BW2 = ((carrier_b_r)/mod_eff)*(1+0.15);% Carrier BW

    % BW2 = 10*log10(carrier_BW2);

    EIRPmax = 39;%Satellite maximum EIRP to receiver location
    Lcontorno=0.5;%Losses due to contour of the link budget
    frequency_downlink = 11114e6;%Higher central frequency for downlink in Ku band
    lambda_link_dw = c/frequency_downlink;%Downlink wavelength
    Diameter_sat = 1; %Antenna Diameter [m]
    min_beam_width_sat =70*lambda_link_dw/Diameter_sat;
    G_down = efficiency*48360/min_beam_width_sat^2;
    G_down = 10*log10(G_down);
    p_HPA_down = 20; %Power of a single HPA multiplied by the n of channels
    p_t_down = 20*log10(p_HPA_down)-l_mc-l_ftx; %Total power for a transponder
    p_tx_down = p_t_down + 10*log10(12); %Total power tx by the antenna (p_t x the n of transponders)
    % EIRPsat = EIRPmax+Lcontorno;%Satellite EIRP
    % percentage_carrier=(Transp_BW/carrier_BW2)*100;
    % EIRP_assigned_to_carrier=10*log10(10^(EIRPmax/10)*(percentage_carrier/100))+OBO;

    Gvsat = 10*log10(efficiency*((pi*Diameter_sat*frequency_downlink)/c)^2);
    Tvsat = 290; %Vsat temperature
    LNB_noise_figure = 1.3;
    T_LNB = 290*((10^(LNB_noise_figure/10))-1);
    GTvsat = Gvsat - 10*log10(Tvsat + T_LNB); %G/T user station
    pointing_loss = 0.3; %[dB]
    gases_absortion = 0.3;
    EIRP_down = G_down + p_tx_down; %Gateway EIRP

    for n = 1:length(t)   
       Ld_pathloss(n) = 20*log10((4*pi*distance_down(n)*frequency_downlink)/c);  %Path Loss
       Ld_rain(n) = (4.21e-5 * (frequency_downlink/10^9)^2.42)*R^(1.41*(frequency_downlink/10^9)^(-0.0779))*2/sind(best_elevation(n)); %I assume that the clouds are 2km high
       CN_downlink(n) = EIRP_down + GTvsat - Ld_pathloss(n) - k - 10*log10(carrier_s_r) - pointing_loss - gases_absortion;
       CN_downlink_rain(n) = CN_downlink(n) - Ld_rain(n);
    end

    %Interferences
    carr_to_inter_noise = 23;

    add_deg = 0.5;
    imp_margin = 1.2;

    required_EsNo=6.55;

    %CN_TOT
    for m = 1:length(t)
    %total link
    CNtot(m) = -20*log10(1/(10^(CN_uplink(m)/20))+1/(10^(CN_downlink(m)/20)+1/(10^(carr_to_inter_noise/20))));
    available_EsNo(m) = CNtot(m) - add_deg - imp_margin;
    link_margin(m) = available_EsNo(m) - required_EsNo;
    
    CNtot_rain(m) = -20*log10(1/(10^(CN_uplink_rain(m)/20))+1/(10^(CN_downlink_rain(m)/20)+1/(10^(carr_to_inter_noise/20))));
    available_EsNo_rain(m) = CNtot_rain(m) - add_deg - imp_margin;
    link_margin_rain(m) = available_EsNo_rain(m) - required_EsNo;
    
    end

end



