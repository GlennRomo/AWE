%-------------------------------------------------------%
%*MainFile*.m                                           %
%                                                       %
%Kite Simulation of a horizontal path, microns off the  %
%ground. The simulation solves for theta, theta_d, l,   %
%l_d, and tension. The simulation involves lift and drag%
%forces and has one animation for the top view of a     %
%flight path with a simple FBD.                         %
%-------------------------------------------------------%

clear all;
close all;
clc;

global D L l_dd Vinfin m Fp C_LiftDragFile B density Vinf_mag S AR u;
global arccomp k wave_count T SYS_time SYS_theta SYS_theta_d SYS_l SYS_l_d;
global LiftVect DragVect Vrel magg_Vkite VkiteVect theta_in theta_out;

%-----------------------------ODE Solution--------------------------------

i = [1,0,0];  %Basic Vector assignment for unit
j = [0,1,0];  %vectors i,j,k
k = [0,0,1];

%Excel Document that contains the lift and drag
%coefficients at 5degree intervals
C_LiftDragFile = xlsread('NACA0015.xlsx');                                  
                                                                            

t_scale = 0:0.01:1000;      %ODE TimeScale
n_waves = 100;

%System Variables
m = 1.2;    %kg          : Mass of the "Kite" (Point Particle).
B = pi()/2;	%deg         : Angle of Chord Vector w/ Tether Vector.
density = 1.225; %kg/m^3 : Wind Power Data Standard [15C @ 1atm].
b = 1.2;	%m           : Wing Span (Based on SenDesign AirFoil)
c = 0.2;	%m           : Chord Length (Based on SenDesign AirFoil)
S = b*c;	%m^2         : Effective Area of Wing
AR = b^2/S;	%            : Aspect Ratio (dimensionless)
u = 1.7976*10^(-5); %kg/(m*s)   : Dynamic Viscosity
Fp = 0;     %N           : A perendicular force to the kite


%Defined Parameters
Vinfin = 10*i;	%m/s    : Wind Velocity.
Vinf_mag = norm(Vinfin); %     : Magnitude of the wind velocity
l_d_out = 2.07;          %m/s  : Length Velocity of Tether during Reel-Out
l_d_in = -3.75;          %m/s  : Length Velocity of Tether during Reel-In
theta_in = 130*pi/180;   %rad  : Theta-In Trigger Location
theta_out = 260*pi/180;  %rad  : Theta-Out Trigger Location
                                                                            
%Initial Conditions
theta = theta_out - 2*pi;%rad  : Angle Position at t=0
theta_d = 0.8761;        %rad/s: Angle Velocity at t=0
l = 118.07;              %m    : Length of Tether at t=0
l_d = l_d_out;           %m/s  : Length Velocity of Tether at t=0
l_dd = 0;                %m/s^2: Acceleration of Tether Release at t=0
D = 0;                   %N    : Drag Force Initial Definition
L = 0;                   %N    : Lift Forces Initial Definition

%System Arrays for Entire Simulation Period
SYS_time = 0;
SYS_theta = theta;
SYS_theta_d = theta_d;
SYS_l = l;
SYS_l_d = l_d;
[~,sizeSYS] = size(SYS_theta);

%Relative and Absolute Tolerances for ODE45 and Event Function Call
my_RelTol = 1e-10;
my_AbsTol = 1e-10;

%SS Tolerance Criteria
ODE_kill = 0;
ODE_kill_end = 2;
SS_tol = 1e-6;

fprintf('theta_d_init: %g l_init: %d \r WaveCount: ', theta_d, l)

%Initial Definition for End Wave Transient Response
End_Wave_l_Trans = l;
End_Wave_theta_d_Trans = theta_d;


%ODE45 Solution Iterations for Theta and L using Event Detection
    % ODE45 Solver with Time and 4 Arc Components: theta, theta_dot,
    % length, and length_dot. The function calls KITE_SS where states
    % are used to incrementally iterate the listed Arc Components

for wave_count = 1:n_waves
    
    
    fprintf('%d ', wave_count);
    

    
    %1st Event ODE45 Solver
        options = odeset('RelTol',my_RelTol,'AbsTol',my_AbsTol,'Events',@PI_2);
        [t1,arccomp] = ode45(@KITE_SS,t_scale,[SYS_theta(sizeSYS),SYS_theta_d(sizeSYS),SYS_l(sizeSYS),SYS_l_d(sizeSYS)],options);
        
        %Partial Wave Fragments from 1st Event Solver
        t2 = t1(:) + SYS_time(sizeSYS);
        time = transpose(t2);theta = transpose(arccomp(:,1));theta_d = transpose(arccomp(:,2));l_wave = transpose(arccomp(:,3));l_d = transpose(arccomp(:,4));
 
        %Last Wave (Pass) Array for State Components
        if ODE_kill > (ODE_kill_end-1)
            LAST_WAVE_time = transpose(t1);
            LAST_WAVE_theta = theta;
            LAST_WAVE_theta_d = theta_d;
            LAST_WAVE_l = l_wave;
            LAST_WAVE_l_d = l_d;
        end
        
        %Total System Arrays addition of Partial Wave Fragment
        SYS_time = [SYS_time time];
        SYS_theta = [SYS_theta theta];
        SYS_theta_d = [SYS_theta_d theta_d];
        SYS_l = [SYS_l l_wave];
        SYS_l_d = [SYS_l_d l_d];
                         
        [~,sizeSYS] = size(SYS_theta);
        [~,sizeL_D] = size(l_d);
        
        %New Controller Input - Change Length Velocity for Reel-In
        SYS_l_d(sizeSYS) = l_d_in;
        
    %Steady State Tolerance Error Definition for Start Wave
        Start_Wave_l(wave_count) = l_wave(1);
        Start_Wave_theta_d(wave_count) = theta_d(1);

        
    %2nd Event ODE45 Solver
        options = odeset('RelTol',my_RelTol,'AbsTol',my_AbsTol,'Events',@PI3_2);
        [t1,arccomp] = ode45(@KITE_SS,t_scale,[SYS_theta(sizeSYS),SYS_theta_d(sizeSYS),SYS_l(sizeSYS),SYS_l_d(sizeSYS)],options); 
        
        %Partial Wave Fragments from 2nd Event Solver
        t2 = t1(:) + SYS_time(sizeSYS);
        time = transpose(t2);theta = transpose(arccomp(:,1));theta_d = transpose(arccomp(:,2));l_wave = transpose(arccomp(:,3));l_d = transpose(arccomp(:,4));

        %Last Wave (Pass) Array for State Components
        if ODE_kill > (ODE_kill_end-1)
            [~,sizeLAST_WAVE] = size(LAST_WAVE_time);
            t3 = t1(:) + LAST_WAVE_time(sizeLAST_WAVE);
            LAST_WAVE_time = [LAST_WAVE_time transpose(t3)];
            LAST_WAVE_theta = [LAST_WAVE_theta theta];
            LAST_WAVE_theta_d = [LAST_WAVE_theta_d theta_d];
            LAST_WAVE_l = [LAST_WAVE_l l_wave];
            LAST_WAVE_l_d = [LAST_WAVE_l_d l_d];
        end
        
        %Total System Arrays addition of Partial Wave Fragment
        SYS_time = [SYS_time time];
        SYS_theta = [SYS_theta theta];
        SYS_theta_d = [SYS_theta_d theta_d];
        SYS_l = [SYS_l l_wave];
        SYS_l_d = [SYS_l_d l_d];
                
        [~,sizeSYS] = size(SYS_theta);
        [~,sizeL_D] = size(l_d);
        
        %New Controller Input - Change Length Velocity for Reel-Out
        SYS_l_d(sizeSYS) = l_d_out;
        
        %End Wave Transient Response
        End_Wave_l_Trans = [End_Wave_l_Trans  l_wave(sizeL_D)];
        End_Wave_theta_d_Trans = [End_Wave_theta_d_Trans  theta_d(sizeL_D)];

     %Steady State Tolerance Error Calculation and End Wave Definition
        End_Wave_l(wave_count) = SYS_l(sizeSYS);
        End_Wave_theta_d(wave_count) = SYS_theta_d(sizeSYS);

        SS_error_l(wave_count) = abs((End_Wave_l(wave_count) - Start_Wave_l(wave_count)) / Start_Wave_l(wave_count));
        SS_error_theta_d(wave_count) = abs((End_Wave_theta_d(wave_count) - Start_Wave_theta_d(wave_count)) / Start_Wave_theta_d(wave_count));
        Total_error(wave_count) = SS_error_l(wave_count) + SS_error_theta_d(wave_count);

        %Steady State Convergence Function Kill Adder
        if Total_error(wave_count) <= SS_tol
            ODE_kill = ODE_kill+1;
        elseif ODE_kill == 1
            ODE_kill = 0;
        end
        
        %Steady State Convergence Function Kill
        if ODE_kill > ODE_kill_end
            break
        elseif ODE_kill > 1
        	ODE_kill = ODE_kill+1;
        end

end

%Determine Minimum Tension of the System
Min_Tension = Tension(1);

%Net Average Cycle Power
LAST_WAVE_Power_avg = Power(1);

%Animation
b = aniMATION(1);
