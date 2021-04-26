function [A, B, d] = rk4_kinematic_curvilinear(x, u, kappa, dt)
%RK4_KINEMATIC_CURVILINEAR Linearises the dynamics of the kinematic
%bicycle model using a curvilinear coordinate frame at a given setpoint
%using Runge-Kutta 4th order scheme
%   INPUTS:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);   
%       dt - Time step
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
    
    rk_factor = [1 2 2 1];
    for i = 1:N_steps
        x_k = x(:, i);
        u_k = u(:, i);
        for j = 1:4
            % Populate matrices  
            A_k = A_curv_kin(x_k, u_k, kappa);
            A(:, :, i) = A(:, :, i) + A_k*rk_factor(j);

            B_k = B_curv_kin(x_k, u_k, kappa);   
            B(:, :, i) = B(:, :, i) + B_k*rk_factor(j);

            f = f_curv_kin(x_k, u_k, kappa);
            d(:, i) = d(:, i) + (f - A_k*x_k - B_k*u_k) * rk_factor(j);

            x_k = x(:, i) + f*dt / rk_factor(mod(j, 4)+1);
        end
    end
    
    A = A / 6;
    B = B / 6;
    d = d / 6;
    
end

