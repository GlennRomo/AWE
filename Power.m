%--------------------------------------------------%
% Power.m                                          %
%                                                  %
% The net average cycle power for a single wave is %
% wanted to know what the system will produce once %
% reaching steady state, periodic motion. This     %
% function will run only if the system convergences%
%--------------------------------------------------%


function [ output ] = Power( input )

  global LAST_WAVE_T LAST_WAVE_l_d LAST_WAVE_time LAST_WAVE_P;

    if ODE_kill > ODE_kill_end
        
        [~,sizeLAST_WAVE] = size(LAST_WAVE_time);

        LAST_WAVE_P = times(LAST_WAVE_T,LAST_WAVE_l_d);
        LAST_WAVE_Energy_total = trapz(LAST_WAVE_time,LAST_WAVE_P);
        LAST_WAVE_Power_avg = LAST_WAVE_Energy_total / LAST_WAVE_time(sizeLAST_WAVE);

        fprintf('Last Wave Energy Average: %.2f W\r', LAST_WAVE_Power_avg);

    end

    output = LAST_WAVE_Power_avg;

end

