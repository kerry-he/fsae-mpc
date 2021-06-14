function [A, B, d] = euler_kinematic_curvilinear(x, u, kappa, ~)
%EULER_KINEMATIC_CURVILINEAR Linearises the dynamics of the kinematic
%bicycle model using a curvilinear coordinate frame at a given setpoint
%using first order Euler integration
%   INPUTS:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);   
%       kappa_d - Anonymous function of spline curvature derivative
%   OUTPUTS:
%       Parameters of the linearised differential equation: 
%           dx/dt = Ax + Bu + d:

    % Define constant sizes
    [N_x, N_steps] = size(x);
    [N_u, ~] = size(u); 

    % Preallocate
    A = zeros(N_x, N_x, N_steps);
    B = zeros(N_x, N_u, N_steps);
    d = zeros(N_x, N_steps);
    
    for i = 1:N_steps
        % Populate matrices
        A(:, :, i) = A_curv_dyn(x(:, i), u(:, i), kappa);
        B(:, :, i) = B_curv_dyn(x(:, i), u(:, i), kappa); 
        f = f_curv_dyn(x(:, i), u(:, i), kappa);
        d(:, i) = f - A(:, :, i)*x(:, i) - B(:, :, i)*u(:, i);
    end
    
end

