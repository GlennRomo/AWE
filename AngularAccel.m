%------------------------------------------------------%
% AngularAccel.m                                       %
%                                                      %
% For-Loop desinged to calculate tension and angular   %
% corresponding to values from the ODE45 found in the  %
% Arc Components array. Tension isn't calculated in the%
% state function because it isn't based on a derivative%
% of itself and is only dependent upon the Arc         %
% Components. Angular Acceleration is not saved during %
% the ODE solver and is therefore recalculated for     %
% incremental time-steps in this function.             %
%------------------------------------------------------%

function [ output ] = AngularAccel( input )

 global SYS_time SYS_theta SYS_theta_d SYS_l SYS_l_d T S density Vinf_mag;
 global l_dd m c u Re LAST_WAVE_T SYS_l_dd;

    %Matrix Size Declarations
    [~,t_size] = size(SYS_time);% : Size of Counter used in For-Loop
    T = zeros(1,t_size);        % : Declare Size of Tension Array
    SYS_l_dd = zeros(1,t_size); % : Declare Size of AngAccel Array
    Re = zeros(1,t_size);       % : Declare Size of Reynold's # Array

    [f,g] = size(SYS_time);
    ang_attack = zeros(1,f);    % : Declare Size of Angle of Attack Array
    gamma = zeros(1,f);         % : Declare Size of Gamma Array

    %Determine Last Wave Dynamics if Simulation Converged
    if ODE_kill > ODE_kill_end
        [~,sizeLAST_WAVE] = size(LAST_WAVE_time);
        LAST_WAVE_count = 1;
    end
    
    %Iterative For-Loop for Size of System Arrays to calculate Tension
    %and Angular Acceleration
    for ind=1:t_size

        %Initial Conditions
        % Define theta, theta-dot, length, length-dot in For-Loop as
        % the respective value from ODE45 Calculations
        theta = SYS_theta(1,ind);
        theta_d = SYS_theta_d(1,ind);
        l = SYS_l(1,ind);
        l_d = SYS_l_d(1,ind);

        %Initial Conditions
            % Define theta, theta-dot, length, length-dot in For-Loop as
            % the respective value from ODE45 Calculations
        [output] = angleOfAttack(theta,theta_d,l,l_d);

        % Angle of Attack Calculations are performed a second time to find
        % values that correspond with theta, theta_d, length, and length_d
        % since the angle isn't stored in an array when evaluated in the
        % variable timestep ODE45 function
        ang_attack(ind) = output(1);
        Vrel(1,ind) = output(2);
        Vrel(2,ind) = output(3);

        DragVect(1,ind) = output(4);
        DragVect(2,ind) = output(5);
        LiftVect(1,ind) = output(6);
        LiftVect(2,ind) = output(7);

        %Coefficient of Lift and Drag Calculation
        [c_DL] = CoeffLiftDrag(ang_attack(ind));

        L(1,ind) = 0.5*S*density*(theta_d*l)^2*c_DL(2);   %N
        D(1,ind) = 0.5*S*density*(theta_d*l)^2*c_DL(1);   %N

        %Kite Velocity Calculation - Necessary for the Animation
        magg_Vkite(ind) = norm([output(8),output(9),output(10)]);   %Component Vectors of Vkite
        VkiteVect(ind,:) = [output(8),output(9),output(10)]/magg_Vkite(ind);

        %Angular Acceleration System Array
        SYS_l_dd(1,ind) = (l*m*theta_d.^2*cos(theta).^2 - SYS_Tension(ind)*sin(theta).^2 - (D(1,ind)*l_d*cos(theta).^2)./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2) - SYS_Tension(ind)*cos(theta).^2 - (D(1,ind)*l_d*sin(theta).^2)./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2) + l*m*theta_d.^2*sin(theta).^2 + (D(1,ind)*Vinf_mag*cos(theta))./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2) + (L(1,ind)*Vinf_mag*sin(theta))./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2) + (L(1,ind)*l*theta_d*cos(theta).^2)./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2) + (L(1,ind)*l*theta_d*sin(theta).^2)./(abs(l_d*sin(theta) + l*theta_d*cos(theta)).^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta)).^2).^(1./2))./(m*(cos(theta).^2 + sin(theta).^2));
        l_dd = SYS_l_dd(1,ind);

        %Tension and Reynolds Number Calculation
        T(1,ind) = (l*m*theta_d^2*cos(theta)^2 - l_dd*m*sin(theta)^2 - (D(1,ind)*l_d*cos(theta)^2)/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) - l_dd*m*cos(theta)^2 - (D(1,ind)*l_d*sin(theta)^2)/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + l*m*theta_d^2*sin(theta)^2 + (D(1,ind)*Vinf_mag*cos(theta))/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + (L(1,ind)*Vinf_mag*sin(theta))/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + (L(1,ind)*l*theta_d*cos(theta)^2)/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + (L(1,ind)*l*theta_d*sin(theta)^2)/(abs(l_d*sin(theta) + l*theta_d*cos(theta))^2 + abs(Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2))/(cos(theta)^2 + sin(theta)^2);
        Re(1,ind) = density*(theta_d*l)*c/u;        

        %If Convergence was Reached
        if ODE_kill > ODE_kill_end
            if ind >= (t_size - sizeLAST_WAVE + 1)
                LAST_WAVE_T(LAST_WAVE_count) = T(1,ind);
                LAST_WAVE_count = LAST_WAVE_count + 1;
            end
        end

    end

    %Minimum and Maximum Tensions
    MAX_Tension = max(T);
    MIN_Tension = min(T);

    fprintf('\nTension - Max: %.2f N  Min: %.2f N\r', MAX_Tension, MIN_Tension);
    
    output = MIN_Tension;

end

