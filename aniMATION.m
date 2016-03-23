%---------------------------------------------------------%
% aniMATION.m                                             %
%                                                         %
% Provides a visual representation of the flightpath and  %
% and orientation of the airfoil while in flight. This    %
% version depicts if the tether is in positive or         %
% negative tension by the color of the thin rod. It will  %
% also provide a FBD with force vectors in ratio to their %
% maximum values.                                         %
%---------------------------------------------------------%

function b = aniMATION(b)

  global SYS_time SYS_theta SYS_theta_d SYS_l SYS_l_d L D LiftVect DragVect;
  global Vrel T magg_Vkite VkiteVect;

    [~,sizeSYS] = size(SYS_time);

    %Initital Defining of Lines
    
        %Open a Figure Window to place initial construct of animation
        verticalFig = figure('Name','2D Tether Animation');                         
        %Line for initial tether length and initial theta
        tether_line = line('xdata',[0 SYS_l(1,1)*cos(SYS_theta(1,1))],'ydata',[0 SYS_l(1,1)*sin(SYS_theta(1,1))]);      
        max_l = max(SYS_l);
        %Tether Figure Margin
        xlim([-max_l-max_l*.2,max_l+max_l*.2]);
        ylim([-max_l-max_l*.2,max_l+max_l*.2]);
        axis equal
        axis manual

    %Unitization of Lift, Drag, and Velocity Lines
        M1 = max(D);
        M2 = max(L);

        if M1 < M2
            maxLD = M2;
        else
            maxLD = M1;
        end

        maxVEL = max(magg_Vkite);

    %Kite Figure
        kiteFig = figure('Name','Kite FBD');
        magg_vrel = norm(Vrel(1,:));
        lift_line = line('xdata',[0 L(1,1)*LiftVect(1,1)/maxLD],'ydata',[0 L(1,1)*LiftVect(1,1)/maxLD]);
        drag_line = line('xdata',[0,D(1,1)*DragVect(1,1)/maxLD],'ydata',[0,D(1,1)*DragVect(1,1)/maxLD]);
        vrel_line = line('xdata',[0,Vrel(1,1)/magg_vrel],'ydata',[0,Vrel(1,2)/magg_vrel]);
        vkite_line = line('xdata',[0,VkiteVect(1,1)/maxVEL],'ydata',[0,VkiteVect(1,2)/maxVEL]);

    %XY-axis limits for Airfoil Figure Window
        xlim([-1.1,1.1]);  
        ylim([-1.1,1.1]);
        axis equal
        axis manual

    %Airfoil Outline Points
        x_foil = [0.00158000000000000;0.0100800000000000;0.0181000000000000;...
            0.0327900000000000;0.0458000000000000;0.0570400000000000;...
            0.0661700000000000;0.0725400000000000;0.0750200000000000;...
            0.0742700000000000;0.0717200000000000;0.0668200000000000;...
            0.0585300000000000;0.0525000000000000;0.0444300000000000;...
            0.0326800000000000;0.0236700000000000;0.00100000000000000;...
            -0.0236700000000000;-0.0326800000000000;-0.0444300000000000;...
            -0.0525000000000000;-0.0585300000000000;-0.0668200000000000;...
            -0.0717200000000000;-0.0742700000000000;-0.0750200000000000;...
            -0.0725400000000000;-0.0661700000000000;-0.0570400000000000;...
            -0.0458000000000000;-0.0327900000000000;-0.0181000000000000;...
            -0.0100800000000000;-0.00158000000000000];
        y_foil = [-0.750000000000000;-0.700000000000000;-0.650000000000000;...
            -0.550000000000000;-0.450000000000000;-0.350000000000000;...
            -0.251000000000000;-0.150000000000000;-0.0500000000000000;...
            0;0.0500000000000000;0.100000000000000;0.150000000000000;...
            0.175000000000000;0.200000000000000;0.225000000000000;...
            0.237500000000000;0.250000000000000;0.237500000000000;...
            0.225000000000000;0.200000000000000;0.175000000000000;...
            0.150000000000000;0.100000000000000;0.0500000000000000;0;...
            -0.0500000000000000;-0.150000000000000;-0.251000000000000;...
            -0.350000000000000;-0.450000000000000;-0.550000000000000;...
            -0.650000000000000;-0.700000000000000;-0.750000000000000];
        
        [q,r] = size(x_foil);

    %Airfoil Shape Array
        for ind = 1:q-1
            airfoil_parts(ind) = line('xdata',[0,0],'ydata',[0,0]);
        end

    %Draw Dynamics and FBD for the entire simulation
        for ind=1:sizeSYS

            %Initial Parameters
            force_scale = 50;
            theta = SYS_theta(1,ind);
            theta_d = SYS_theta_d(1,ind);
            l = SYS_l(1,ind);
            l_d = SYS_l_d(1,ind);
            Lift = L(1,ind)/maxLD;
            Drag = D(1,ind)/maxLD;
            Vel = magg_Vkite(ind)/maxVEL;

            %Rotate and Draw Airfoil for each iteration
            for indd = 1:q
                [foilCoord1(:,indd)] = [cos(theta),-sin(theta);sin(theta),cos(theta)]*[x_foil(indd,1);y_foil(indd,1)];       
            end

            foilCoord = foilCoord1';

            for indd = 1:q-1
                set(airfoil_parts(indd),'xdata',[foilCoord(indd,1),foilCoord(indd+1,1)],'ydata',[foilCoord(indd,2),foilCoord(indd+1,2)]);
            end

            % 2-D Line Drawing for the Simulation
                %Tether Figure (verticalFig)
                if T(1,ind)>0
                    set(tether_line,'xdata',[0,l*cos(theta)],'ydata', [0,l*sin(theta)]);
                else
                    set(tether_line,'xdata',[0,l*cos(theta)],'ydata', [0,l*sin(theta)],'Color','red');
                end

                %Airfoil Lift and Drag Components Animation (kiteFig)
                mag_vrel = norm(Vrel(:,ind));
                set(vrel_line,'xdata',[0,Vrel(1,ind)/mag_vrel],'ydata',[0,Vrel(2,ind)/mag_vrel],'Color','yellow');
                set(lift_line,'xdata',[0,Lift*LiftVect(1,ind)],'ydata',[0,Lift*LiftVect(2,ind)],'Color','green');
                set(drag_line,'xdata',[0,Drag*DragVect(1,ind)],'ydata',[0,Drag*DragVect(2,ind)],'Color','red');
                set(vkite_line,'xdata',[0,Vel*VkiteVect(ind,1)],'ydata',[0,Vel*VkiteVect(ind,2)],'Color','blue');

            drawnow
            grid on

end