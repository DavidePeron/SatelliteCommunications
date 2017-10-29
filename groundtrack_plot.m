function groundtrack_plot(earthmap, t, r, re, El, long, lat, color, lat_gs, long_gs)

    n_sat = length(long);
    h1 = cell(n_sat);
    for i=0:n_sat-1
        h1{i+1} = plot(0,0,'o','Markersize',1);
    end
    
    figure (2);
    set(gcf,'Menubar','none','Name','Earth Track', ... 
        'NumberTitle','off','Position',[70,30,1000,500]); 
    hold on
    image([180 -180],[90 -90],earthmap,'CDataMapping', 'scaled');
    axis equal
    axis ([-180 180 -90 90]);
    set(gca, 'XDir', 'reverse');

    %% Ground Stations plot
    for i=1:length(long_gs)    
        plot (long_gs(i),lat_gs(i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 5);
%         set (gs(i), 'Visible', 'on');
    end
    
    for i=1:length(t)
        for j = 1:n_sat
            %% Nadir angle estimation
            rho = asin(re/r{j}(i));
            beta = asin(cos(deg2rad(El))*sin(rho));
            nadir = beta;

            %% Lambda and arch estimation 
            lambda= 90 - rad2deg(nadir) - El; %central angle in degrees

            plot (long{j}(i),lat{j}(i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',color(j),'MarkerSize', 2);
            if (i~=1 && abs(long{j}(i-1)-long{j}(i))<100)
                line([long{j}(i-1) long{j}(i)],[lat{j}(i-1) lat{j}(i)],'Color', color(j), 'LineWidth', 2);
            end
            delete(h1{j});
            h1{j} = plot(long{j}(i),lat{j}(i),'o','MarkerEdgeColor',color(j), 'Markersize', lambda*2);
        end

        pause (0.01);
    end

end