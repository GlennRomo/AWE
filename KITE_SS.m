%-------------------------------------------------------%
%KITE_SS.m                                              %
%                                                       %
%This file takes an input from *MainFile*.m to          %
%produce a numerical solution for four variables: theta,%
%theta_d, length, and length_d through the use of a     %
%State Space Model.                                     %
%-------------------------------------------------------%



function [myoutput,T] = KITE_SS(t,state)
    
    global D L l_dd Vinfin m Fp ang_attack B density Vinf_mag S;
    

    
    %Initial State Declarations
    theta = state(1);                                                       
    theta_d = state(2);                                                     
    l = state(3);                                                           
    l_d = state(4);  
    
    % The angle of attack function is called sending the state-variables to
    % find the angle of attack in the current state. Angle must also be
    % sent in radians.                                      
    
    %Read the Excel Sheet for Coeff of Lift and Drag
    output = angleOfAttack(theta,theta_d,l,l_d);
    ang_attack = output(1);                                                                        
    [c_DL] = CoeffLiftDrag(ang_attack);
    
    L = 0.5*S*density*(theta_d*l)^2*c_DL(2);
    D = 0.5*S*density*(theta_d*l)^2*c_DL(1);
        
    % State Variables that depend on derivatives. These are taken from the
    % KiteSimul file as initial conditions and then incrementally found by
    % plugging the initials values into the equations for derivatives and
    % original with respect to the ODE45 solver. This solver will find the
    % value evaluated at variable stepsizes depending on activity of the
    % function.
    
    %State Changes Based on EOM
    state_dot(1) = theta_d;
    state_dot(2) = -((D*Vinf_mag*sin(theta))/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) - Fp*sin(theta)^2 - (L*Vinf_mag*cos(theta))/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) - Fp*cos(theta)^2 + (L*l_d*cos(theta)^2)/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + (L*l_d*sin(theta)^2)/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + 2*l_d*m*theta_d*cos(theta)^2 + 2*l_d*m*theta_d*sin(theta)^2 + (D*l*theta_d*cos(theta)^2)/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2) + (D*l*theta_d*sin(theta)^2)/((l_d*sin(theta) + l*theta_d*cos(theta))^2 + (Vinf_mag - l_d*cos(theta) + l*theta_d*sin(theta))^2)^(1/2))/(l*m*(cos(theta)^2 + sin(theta)^2));
    state_dot(3) = l_d;
    state_dot(4) = 0;
    
    myoutput = [state_dot(1);state_dot(2);state_dot(3);state_dot(4)];


    
    