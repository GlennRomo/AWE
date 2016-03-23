%-------------------------------------------------------%
%angleOfAttack.m                                        %
%                                                       %
%An input is taken from either *MainFile.m* or Kite_SS.m%
%to find the angle of attack based on theta, theta_d, l,%
%and l_d. The angle is calculated by subtracting gamma  %
%from our angle Beta the kite makes with the tether.    %
%-------------------------------------------------------%

function [output] = angleOfAttack(theta,theta_d,l,l_d)

global Vinfin B k;

%Er and Etheta components in terms of i,j,k
e_r = [cos(theta),sin(theta),0];
e_theta = [-sin(theta),cos(theta),0];

%Velocity and Acceleration Components
Vi = Vinfin;
Vkite = theta_d*l*e_theta+l_d*e_r;
Vrel = Vi-Vkite;

%Solving for Angle of Attack
k_vrelteth = (-Vrel(1))*(-e_r(2)) - (-Vrel(2))*(-e_r(1));
Vrel_unit = norm(-Vrel);
Teth_unit = norm(-e_r);
gamma = acos(dot(-Vrel,-e_r)/(Vrel_unit*Teth_unit));    

%Unit Vectors of Lift and Drag
lambda_D = Vrel/Vrel_unit;
lambda_L = cross(k,lambda_D);

%Gamma Correction
if k_vrelteth < 0
    gamma = -gamma;
end

%Solver for Gamma
ang_attack = B - gamma;
ang_attack = wrapTo2Pi(ang_attack);

if ang_attack > pi
    ang_attack = ang_attack - 2*pi;
end    

output = [ang_attack,Vrel(1),Vrel(2),lambda_D(1),lambda_D(2),lambda_L(1),lambda_L(2),Vkite(1),Vkite(2),Vkite(3)];
