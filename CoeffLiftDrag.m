%-------------------------------------------------------%
%CoeffLiftDrag.m                                        %
%                                                       %
%An input for angle of attack is used to interpolate the%
%excel sheet opened in *MainFile*.m to find the         %
%coefficients of lift and drag. The excel sheet has     %
%corrected lift and drag coefficients based on the      %
%airfoil NACA0015.                                      %
%-------------------------------------------------------%


function [myoutput] = CoeffLiftDrag(ang_attck)      
        %Accept Angle of Attack in Radians

global C_LiftDragFile

%Interpolation for Coefficient of Drag and Lift
c_drag = interp1q(C_LiftDragFile(:,2),C_LiftDragFile(:,6),ang_attck);
c_lift = interp1q(C_LiftDragFile(:,2),C_LiftDragFile(:,5),ang_attck);

myoutput = [c_drag,c_lift];