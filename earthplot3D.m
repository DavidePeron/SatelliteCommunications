function earthplot3D(earthmap, t, e, a, ECI, color)

    load('earth_constants.mat');
    n_sat = length(ECI);
    
    set(gcf,'Menubar','default','Name','Orbit Visualization', ... 
        'NumberTitle','off','Position',[70,10,750,750]); 
    lim=(1+e)*a;%Setting the limits of the graph
    clf
    axis([-lim, lim, -lim, lim, -lim, lim])	
    view(150,15) 
    axis equal
    shg
    hold on
    grid on
    title('Orbital Visualization');

    %Plotting the Earth
    [xx yy zz]=ellipsoid (0,0,0,a_WGS84, a_WGS84, b);
    earthmap2 = flip(earthmap,1);
    earthmap2 = flip(earthmap2,2);
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';
    pro.Cdata = earthmap2;
    earth= surface(xx,yy,zz,pro);

    rotate (earth, [0 0 1], 0);
    Xaxis= line([-1e7 1e7],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
    Yaxis= line([0 0],[-1e7 1e7],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
    rotate (Xaxis, [0 0 1], 0);
    rotate (Yaxis, [0 0 1], 0);
    % Sun=light('Position',[1 0 0],'Style','infinite');

    %Plotting the ECI Axes
    line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
    line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
    line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')

    sphere_position = cell(n_sat);
    position = cell(n_sat);

    for i=1:length(t)
        rotate (earth, [0 0 1], (360/length(t)), [0,0,0]);
        rotate (Xaxis, [0 0 1], (360/length(t)), [0,0,0]);
        rotate (Yaxis, [0 0 1], (360/length(t)), [0,0,0]);

        %Drawing the red sphere
        for j=1:n_sat
            %Drawing the red sphere
            sphere_position{j}(i)=plot3 (ECI{j}(1,i), ECI{j}(2,i), ECI{j}(3,i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor',color(j),'MarkerSize', 6);
            position{j}(i)=line([0 ECI{j}(1,i)],[0 ECI{j}(2,i)], [0 ECI{j}(3,i)],'Color', color(j), 'LineWidth', 2);

            if(i~=1)
                set (sphere_position{j}(i-1), 'Visible', 'off');
                set (position{j}(i-1), 'Visible', 'off');
            end

            if (i~=1 && i<=length(t))
                line([ECI{j}(1,i-1) ECI{j}(1,i)],[ECI{j}(2,i-1) ECI{j}(2,i)], [ECI{j}(3,i-1) ECI{j}(3,i)],'Color', 'black', 'LineWidth', 1.5);
            end
        end

        pause (0.01);
    end

end