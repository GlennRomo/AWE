%-------------------------------------------------------%
%Kite_Simul_0216.m                                      %
%                                                       %
%Kite Simulation of a horizontal path, microns off the  %
%ground. The simulation solves for theta, theta_d, l,   %
%l_d, and tension. The simulation involves lift and drag%
%forces and has one animation for the top view of a     %
%flight path.                                           %
%-------------------------------------------------------%

%***------------------------------------------------------------------***%
%***------------Insert Initial Conditions and Parameters--------------***%
%***------------------------------------------------------------------***%

x_star = l_init;

SS_tol = 1e-7; %Tolerance level set for finding a Zero Crossing
del = 0.0001; %Delta Change for Equation Slope

x_star_count = 0; %Completed Passes of Root Find
End_Wave_l = 100000; %Initital Difference between Beginning and End of Wave
Start_Wave_l = 0; %Initial guess for the Root Find
Total_error = 10; %Error Initial Value
        
%Iterative Process Based on the Difference between the beginning of the
%wave and the end of the wave for Angular Velocity or Tether Length are
%both less than the SteadyState Tolerance Level
while abs(End_Wave_l - Start_Wave_l) > SS_tol
    x_star_count = x_star_count+1;

    %Root Find Initial Conditions
    SYS_time = 0;
    SYS_theta = theta;
    SYS_theta_d = theta_d_init;
    SYS_l = x_star(x_star_count,1);
    SYS_l_d = l_d;

    %Initial Solver for Baseline of Root Find
    baseODE = ODEventFuncDIFF();
    Start_Wave_l = baseODE(1); Start_Wave_theta_d = baseODE(2); l_d_out_base = baseODE(3);
    End_Wave_l = baseODE(4); End_Wave_theta_d = baseODE(5); l_d_in_base = baseODE(6);
    Total_error = baseODE(7); theta_in_base = theta_in; theta_out_base = theta_out;

    if abs(End_Wave_l - Start_Wave_l) > SS_tol

        %Define Upper Bound
        del_l = Start_Wave_l*del;
        Start_Wave_l_upper = Start_Wave_l + del_l;

        %Del L Integration to determine slope from Baseline and Upper Bound
        SYS_time = 0;
        SYS_theta = theta;
        SYS_theta_d = Start_Wave_theta_d;
        SYS_l = Start_Wave_l_upper;
        SYS_l_d = l_d;

        del_l_ODE = ODEventFuncDIFF();
        Start_Wave_l_del_l = del_l_ODE(1); Start_Wave_theta_d_del_l = del_l_ODE(2); l_d_out_del_l = del_l_ODE(3);
        End_Wave_l_del_l = del_l_ODE(4); End_Wave_theta_d_del_l = del_l_ODE(5); l_d_in_del_l = del_l_ODE(6);
        Total_error_del_l = del_l_ODE(7); theta_in_del_l = theta_in; theta_out_del_l = theta_out;

        %Differences between Baseline and Upper Bound start and end lengths
        Func_l = (End_Wave_l - Start_Wave_l);
        Func_l_del_l = (End_Wave_l_del_l - Start_Wave_l_del_l);

        %Alpha Coefficient used to calculate the size of the step
        l_der = (Func_l_del_l - Func_l)/del_l;

        %Stepped Value location fed to next loop iteration
        stepped_value = (1/l_der)*-(End_Wave_l-Start_Wave_l) + Start_Wave_l;
        x_star(x_star_count+1,1) = stepped_value(1);

    end

end

