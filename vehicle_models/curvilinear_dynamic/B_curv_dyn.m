function B = B_curv_kin(~, ~, ~)
%B_CURV_KIN Calculates the Jacobian of f_curv_kin with respect to the
%control variables
%   INPUT:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);  
%   OUTPUT:
%       B - Jacobian of vehicle model with respect to control variables

    K_steer = 5.0;

    B = [0   0;
         0   0;
         0   0;
         1   0;
         0   0;
         0   0;
         0   K_steer];

end

