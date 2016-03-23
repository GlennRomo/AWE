%-------------------------------------------------------%
%KITE_ThetaSol.m                                        %
%                                                       %
%This funciton provides function solutions for angular  %
%acceleration and tension. This version has a constant  %
%force applied.                                         %
%-------------------------------------------------------%

i = [1,0,0];
j = [0,1,0];
k = [0,0,1];

%Symbolic Variables
syms T Vinfin D L m Fp B ang_attck
syms l l_d l_dd
syms theta theta_d theta_dd

%Er and Etheta components in terms of i,j,k
e_r = [cos(theta),sin(theta),0];
e_theta = [-sin(theta),cos(theta),0];

%Velocity and Acceleration Components
Vi = Vinfin*i;
Vkite = theta_d*l*e_theta+l_d*e_r;
Vrel = Vi-Vkite;

Akite = (theta_dd*l+2*theta_d*l_d)*e_theta+(l_dd-theta_d^2*l)*e_r;

%Solving for Angle of Attack
k_vrelteth = cross(Vrel,-e_r);
Vrel_unit = norm(Vrel);
Teth_unit = norm(-e_r);
gamma = asin(norm(k_vrelteth)/(Vrel_unit*Teth_unit));

%Unit Vectors of Lift and Drag
lambda_D = Vrel/Vrel_unit;
lambda_L = cross(k,lambda_D);

%Solver for Tension and Theta_dd
eq1 = Fp*e_theta-T*e_r+D*lambda_D+L*lambda_L-m*Akite;
[T,theta_dd] = solve(eq1(1),eq1(2),'T','theta_dd');

%Solver for Gamma
eq2 = B - gamma - ang_attck;
[ang_attck] = solve(eq2(1),'ang_attck');
    
    
    