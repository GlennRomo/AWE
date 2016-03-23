%------------------------------------------------------------%
% 2DRootFind_theta_d_l.m                                     %
%                                                            %
% Newtonian-Rhapsonian Multi-Dimensional Root Finding Method %
% used to find the error plane crossing for the variables    %
% angular velocity and tether length at the end of a cycle.  %
%------------------------------------------------------------%

%***------------------------------------------------------------------***%
%***------------Insert Initial Conditions and Parameters--------------***%
%***------------------------------------------------------------------***%

%Initital Difference between Beginning and End of Wave
Func_theta_d = 100000;
Func_l = 100000;

SS_tol = 1e-5; %Tolerance level set for finding a Zero Crossing
del = 0.0001;  %Delta Change for Equation Slope
x_star = [theta_d, l]; %Initial guess for the Root Find
x_star_count = 0; %Completed Passes of Root Find

%Iterative Process Based on the Difference between the beginning of the
%wave and the end of the wave for Angular Velocity and Tether Length are
%both less than the SteadyState Tolerance Level
while abs(Func_theta_d) > SS_tol || abs(Func_l) > SS_tol
    x_star_count = x_star_count+1;

    %Root Find Initial Conditions
    SYS_time = 0;
    SYS_theta = theta;
    SYS_theta_d = x_star(x_star_count,1);
    SYS_l = x_star(x_star_count,2);
    SYS_l_d = l_d;
    
    %Initial Solver for Baseline of Root Find
    baseODE = ODEventFuncDIFF();
    Start_Wave_l = baseODE(1); Start_Wave_theta_d = baseODE(2); l_d_out_base = baseODE(3);
    End_Wave_l = baseODE(4); End_Wave_theta_d = baseODE(5); l_d_in_base = baseODE(6);
    Total_error = baseODE(7); theta_in_base = theta_in; theta_out_base = theta_out;
    
    %Will enter RF if error between variables do not satisfy SS_tolerance
    if abs(Func_theta_d) > SS_tol || abs(Func_l) > SS_tol
        
        %Define Upper Bound for Angular Velocity and Tether Length
        del_l = Start_Wave_l*del;
        del_theta_d = Start_Wave_theta_d*del;
        
        Start_Wave_l_upper = Start_Wave_l + del_l;
        Start_Wave_theta_d_upper = Start_Wave_theta_d + del_theta_d;
        
        %Del Theta_d Integration
            %Initial Conditions for Del Theta_D
            SYS_time = 0;
            SYS_theta = theta;
            SYS_theta_d = Start_Wave_theta_d_upper;
            SYS_l = Start_Wave_l;
            SYS_l_d = l_d;
            
            %Theta_D ODE Solver
            del_theta_d_ODE = ODEventFuncDIFF();
            Start_Wave_l_del_t = del_theta_d_ODE(1); Start_Wave_theta_d_del_t = del_theta_d_ODE(2); l_d_out_del_t = del_theta_d_ODE(3);
            End_Wave_l_del_t = del_theta_d_ODE(4); End_Wave_theta_d_del_t = del_theta_d_ODE(5); l_d_in_del_t = del_theta_d_ODE(6);
            Total_error_del_t = del_theta_d_ODE(7); theta_in_del_t = theta_in; theta_out_del_t = theta_out;

        %Del L Integration
            %Initial Conditions for Del L
            SYS_time = 0;
            SYS_theta = theta;
            SYS_theta_d = Start_Wave_theta_d;
            SYS_l = Start_Wave_l_upper;
            SYS_l_d = l_d;

            %L ODE Solver
            del_l_ODE = ODEventFuncDIFF();
            Start_Wave_l_del_l = del_l_ODE(1); Start_Wave_theta_d_del_l = del_l_ODE(2); l_d_out_del_l = del_l_ODE(3);
            End_Wave_l_del_l = del_l_ODE(4); End_Wave_theta_d_del_l = del_l_ODE(5); l_d_in_del_l = del_l_ODE(6);
            Total_error_del_l = del_l_ODE(7); theta_in_del_l = theta_in; theta_out_del_l = theta_out;


        %Jacobian, StepSize, and Updated Root
            %Function Calculations for Differences between Begin and End
            Func_theta_d = (End_Wave_theta_d - Start_Wave_theta_d);
            Func_l = (End_Wave_l - Start_Wave_l);
            Func_theta_d_del_t = (End_Wave_theta_d_del_t - Start_Wave_theta_d_del_t);
            Func_theta_d_del_l = (End_Wave_theta_d_del_l - Start_Wave_theta_d_del_l);
            Func_l_del_t = (End_Wave_l_del_t - Start_Wave_l_del_t);
            Func_l_del_l = (End_Wave_l_del_l - Start_Wave_l_del_l);

            %Newton-Rhapson Method Jacobian
            Jacobian = [(Func_theta_d_del_t - Func_theta_d)/del_theta_d, (Func_theta_d_del_l - Func_theta_d)/del_l;
                        (Func_l_del_t - Func_l)/del_theta_d, (Func_l_del_l - Func_l)/del_l];

            %New Root Step for the Next Iteration
            stepped_value = inv(Jacobian)*[-(End_Wave_theta_d-Start_Wave_theta_d);-(End_Wave_l-Start_Wave_l)] + [Start_Wave_theta_d;Start_Wave_l];            
            x_star(x_star_count+1,1) = stepped_value(1); x_star(x_star_count+1,2) = stepped_value(2);

    end
 
end

  