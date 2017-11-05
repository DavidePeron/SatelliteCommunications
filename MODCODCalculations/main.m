clear all;
close all;

% 52 Mhz Transponders
n_trans = 49;
B_trans = 52e6;

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
    Rsymb_downlink = 50e6/nu;
    B_downlink = Rsymb_downlink*(1+0.15);
    if(B_downlink < B_trans)
        Rsymb_uplink = 5e6/nu;
        B_uplink = Rsymb_uplink*(1+0.15);
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
        Rsymb_downlink = 50e6/nu;
        B_downlink = Rsymb_downlink*(1+0.15);
        if(B_downlink < B_trans)
            Rsymb_uplink = 5e6/nu;
            B_uplink = Rsymb_uplink*(1+0.15);
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
 
Rsymb_downlink = 50e6/nu;
B_downlink = Rsymb_downlink*(1+0.15);
Rsymb_uplink = 5e6/nu;
B_uplink = Rsymb_uplink*(1+0.15);
% End of the first solution

%%