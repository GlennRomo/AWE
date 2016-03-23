%-------------------------------------------------------%
% *MainFile*.m                                          %
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

global D L l_dd Vinfin m Fp C_LiftDragFile B density Vinf_mag S AR u k;
global SYS_time SYS_theta SYS_theta_d SYS_l SYS_l_d LiftVect DragVect Vrel;
global theta_in theta_out arccomp wave_count Tension magg_Vkite VkiteVect;
global SYS_l_dd;

%-----------------------------ODE Solution--------------------------------

i = [1,0,0];  %Basic Vector assignment for unit
j = [0,1,0];  %vectors i,j,k
k = [0,0,1];

%Excel Document that contains the lift and drag
%coefficients at 5degree intervals
C_LiftDragFile = xlsread('NACA0015.xlsx');
                                                                            

t_scale = 0:0.01:1000;  %ODE TimeScale
n_waves = 50;           %Maximum Number of Passes for AWE

%System Variables
m = 1.5;     %kg         : Mass of the "Kite" (Point Particle)
B = pi()/2;  %deg        : Angle of Chord Vector w/ Tether Vector.
BIN = pi/2;  %deg        : Angle of Chord and Tether for Reel-In
BOUT = pi/2; %deg        : Angle of Chord and Tether for Reel-Out
density = 1.225; %kg/m^3 : Wind Power Data Standard [15C @ 1atm].
b = 1.2;     %m          : Wing Span (Based on SenDesign AirFoil)
c = 0.2;     %m          : Chord Length (Based on SenDesign AirFoil)
S = b*c;     %m^2        : Effective Area of Wing
AR = b^2/S;  %           : Aspect Ratio (dimensionless)
u = 1.7976*10^(-5); %kg/(m*s)   : Dynamic Viscosity 
Fp = 0;      %N          : A perendicular force to the kite

%Defined Parameters
Vinfin = 10*i;           %m/s : Wind Velocity.
Vinf_mag = norm(Vinfin); %m/s : Magnitude of the wind velocity.
l_d_out = 0;             %m/s : Length Velocity of Tether during Reel-Out
%l_d_in = -1;            %m/s : Length Velocity of Tether during Reel-In
theta_in = 80*pi/180;    %rad : Theta-In Trigger Location
theta_out = 270*pi/180;  %rad : Theta-Out Trigger Location
TensionOUT = 540;        %N   : Set Force in Tether during Reel-Out
TensionIN = 35;          %N   : Set Force in Tether during Reel-In

%Initial Conditions
theta = theta_out - 2*pi; %rad      : Angle Position at t=0
theta_d = 1;              %rad/s    : Angle Velocity at t=0
l = 50;                   %m        : Length of Tether at t=0
l_d = l_d_out;            %m/s      : Length Velocity of Tether at t=0
l_dd = 0;                 %m/s^2    : Acceleration of Tether Release at t=0
D = 0;                    %N        : Drag Force Initial Definition
L = 0;                    %N        : Lift Force Initial Definition

%System Arrays for Entire Simulation Period                                                                            %             components.
SYS_time = 0;
SYS_theta = theta;
SYS_theta_d = theta_d;
SYS_l = l;
SYS_l_d = l_d;
SYS_Tension = Tension;
[~,sizeSYS] = size(SYS_theta);

%Relative and Absolute Tolerances for ODE45 and Event Function Call
my_RelTol = 1e-10;
my_AbsTol = 1e-10;

%SS Tolerance Criteria
ODE_kill = 0;
ODE_kill_end = 2;
SS_tol = 1e-3;

%Initial Condition Figure Title Creation
figThIN = num2str(theta_in*180/pi);
figThOUT = num2str(theta_out*180/pi);
figThD = num2str(theta_d);
figL = num2str(l);
figBIN = num2str(BIN*180/pi);
figBOUT = num2str(BOUT*180/pi);

SSFigure_TITLE = strcat({'SteadyState__ThIN_'},figThIN,{'_ThOUT_'},figThOUT,{'_ThD_'},figThD,{'_L_'},figL,{'_BIN_'},figBIN,{'_BOUT_'},figBOUT);
SSFigure_TITLE = char(SSFigure_TITLE(1));
fprintf('%s \r WaveCount: ',SSFigure_TITLE)


%ODE45 Solution Iterations for Theta and L using Event Detection
    % ODE45 Solver with Time and 4 Arc Components: theta, theta_dot,
    % length, and length_dot. The function calls KITE_SS where states
    % are used to incrementally iterate the listed Arc Components
    
for wave_count = 1:n_waves
    
    
    fprintf('%d ', wave_count);
    
    Tension = TensionOUT;
    
    %1st Event ODE45 Solver
        options = odeset('RelTol',my_RelTol,'AbsTol',my_AbsTol,'Events',@PI_2);
        [t1,arccomp] = ode45(@KITE_SS,t_scale,[SYS_theta(sizeSYS),SYS_theta_d(sizeSYS),SYS_l(sizeSYS),SYS_l_d(sizeSYS)],options);
        
        %Partial Wave Fragments from 1st Event Solver
        t2 = t1(:) + SYS_time(sizeSYS);
        time = transpose(t2);theta = transpose(arccomp(:,1));theta_d = transpose(arccomp(:,2));l = transpose(arccomp(:,3));l_d = transpose(arccomp(:,4));
        [~,sizeWAVE] = size(time);
        
        %Last Wave (Pass) Array for State Components
        if ODE_kill > (ODE_kill_end-1)
            LAST_WAVE_time = transpose(t1);
            LAST_WAVE_theta = theta;
            LAST_WAVE_theta_d = theta_d;
            LAST_WAVE_l = l;
            LAST_WAVE_l_d = l_d;
            
            LAST_WAVE_Tension = ones(1,sizeWAVE)*Tension;
        end
        
        %Total System Arrays addition of Partial Wave Fragment
        if wave_count == 1
            SYS_time = time;
            SYS_theta = theta;
            SYS_theta_d = theta_d;
            SYS_l = l;
            SYS_l_d = l_d;

            [~,sizeSYS] = size(SYS_theta);

            SYS_Tension = ones(1,sizeWAVE)*Tension;
        end
            
        if wave_count > 1
            SYS_time = [SYS_time time];
            SYS_theta = [SYS_theta theta];
            SYS_theta_d = [SYS_theta_d theta_d];
            SYS_l = [SYS_l l];
            SYS_l_d = [SYS_l_d l_d];

            [~,sizeSYS] = size(SYS_theta);


            Wave_Tension = ones(1,sizeWAVE)*Tension;
            SYS_Tension = [SYS_Tension Wave_Tension];
        end
        
        
    %Steady State Tolerance Error Definition for Start Wave
        Start_Wave_l(wave_count) = l(1);
        Start_Wave_theta_d(wave_count) = theta_d(1);

    %Tension Switch for New Controller Input
    Tension = TensionIN;
        
    %2nd Event ODE45 Solver
        options = odeset('RelTol',my_RelTol,'AbsTol',my_AbsTol,'Events',@PI3_2);
        [t1,arccomp] = ode45(@KITE_SS,t_scale,[SYS_theta(sizeSYS),SYS_theta_d(sizeSYS),SYS_l(sizeSYS),SYS_l_d(sizeSYS)],options); 
        
        %Partial Wave Fragments from 2nd Event Solver
        t2 = t1(:) + SYS_time(sizeSYS);
        time = transpose(t2);theta = transpose(arccomp(:,1));theta_d = transpose(arccomp(:,2));l = transpose(arccomp(:,3));l_d = transpose(arccomp(:,4));
        [~,sizeWAVE] = size(time);
        
        %Last Wave (Pass) Array for State Components
        if ODE_kill > (ODE_kill_end-1)
            [~,sizeLAST_WAVE] = size(LAST_WAVE_time);
            t3 = t1(:) + LAST_WAVE_time(sizeLAST_WAVE);
            LAST_WAVE_time = [LAST_WAVE_time transpose(t3)];
            LAST_WAVE_theta = [LAST_WAVE_theta theta];
            LAST_WAVE_theta_d = [LAST_WAVE_theta_d theta_d];
            LAST_WAVE_l = [LAST_WAVE_l l];
            LAST_WAVE_l_d = [LAST_WAVE_l_d l_d];
            
            Tension_WAVE2 = ones(1,sizeWAVE)*Tension;
            LAST_WAVE_Tension = [LAST_WAVE_Tension Tension_WAVE2];
        end
        
        %Total System Arrays addition of Partial Wave Fragment
        SYS_time = [SYS_time time];
        SYS_theta = [SYS_theta theta];
        SYS_theta_d = [SYS_theta_d theta_d];
        SYS_l = [SYS_l l];
        SYS_l_d = [SYS_l_d l_d];
                
        [~,sizeSYS] = size(SYS_theta);

        Wave_Tension = ones(1,sizeWAVE)*Tension;
        SYS_Tension = [SYS_Tension Wave_Tension];

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
 
Min_Tension = AngularAccel(1);

LAST_WAVE_Power_avg = Power(1);

