function [r_t, phi_t] = get_polar(t, T, e, a, mu, n_rev)
    t_p = 0:T:T*n_rev; % time at perigee pass at each cycle
    eta_0 = 2*pi/T; %Angular velocity of the fictitious satellite [rad/sec]
    r_t = zeros(1,length(t));
    v_t = zeros(1,length(t));
    phi_t = zeros(1,length(t));
    M = zeros(1,length(t));
    E = zeros(1,length(t));
    %Computation of polar coordinates 
    for i = 1:length(t)

        flag = 0;
        %Mean Anomaly
        M(i) = eta_0*(t(i) - t_p(1,floor(t(i)/T)+1));

        % Computation of Eccentric Anomaly with Newton-Raphson Method
        E(i) = pi;

        %Newton-Raphson Method
        while flag == 0
            temp_1 = E(i);
            E(i) = E(i) + (M(i) + e*sin(E(i)) - E(i))/(1 - e*cos(E(i)));
            temp_2 = E(i);
            if abs(temp_1 - temp_2) < 1e-2
                flag = 1;
            end
        end


        r_t(i) = a*(1-e*cos(E(i)));
        v_t(i) = sqrt(mu*(2/r_t(i)-1/a));
        phi_t(i) = acos(1/e*((a*(1-e^2))/r_t(i)-1));

        if(E(i) > pi)
            phi_t(i) = 2*pi - phi_t(i);
        end
    end
end