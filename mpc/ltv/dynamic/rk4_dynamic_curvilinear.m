function [A, B, d] = rk4_dynamic_curvilinear(x, u, kappa, dt)
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
    I = eye(N_x);
    A = zeros(N_x, N_x, N_steps);
    B = zeros(N_x, N_u, N_steps);
    d = zeros(N_x, N_steps);
    
    for i = 1:N_steps
        x_i = x(:, i);
        u_i = u(:, i);
        
        % Calculate RK4 slopes
        k1 = f_curv_dyn(x_i, u_i, kappa);
        k2 = f_curv_dyn(x_i + k1*dt / 2, u_i, kappa);
        k3 = f_curv_dyn(x_i + k2*dt / 2, u_i, kappa);
        k4 = f_curv_dyn(x_i + k3*dt, u_i, kappa);
        
        f = (k1 + 2*k2 + 2*k3 + k4) / 6;
        
        % Calculate RK4 partial derivatives w.r.t. states 
        dfdx1 = A_curv_dyn(x_i, u_i, kappa);
        dfdx2 = A_curv_dyn(x_i + k1*dt / 2, u_i, kappa);
        dfdx3 = A_curv_dyn(x_i + k2*dt / 2, u_i, kappa);
        dfdx4 = A_curv_dyn(x_i + k3*dt, u_i, kappa);
        
        dkdx1 = dfdx1;
        dkdx2 = dfdx2 * (I + dkdx1*dt / 2);
        dkdx3 = dfdx3 * (I + dkdx2*dt / 2);
        dkdx4 = dfdx4 * (I + dkdx3*dt);
        
        % Calculate RK4 partial derivatives w.r.t. controls
        dkdu1 = B_curv_dyn(x_i, u_i, kappa);
        dkdu2 = B_curv_dyn(x_i + k1*dt / 2, u_i, kappa) + dfdx2*dkdu1*dt / 2;
        dkdu3 = B_curv_dyn(x_i + k2*dt / 2, u_i, kappa) + dfdx3*dkdu2*dt / 2;
        dkdu4 = B_curv_dyn(x_i + k3*dt, u_i, kappa) + dfdx4*dkdu3*dt / 2;
        
        % Calculate final matrices
        A(:, :, i) = (dkdx1 + 2*dkdx2 + 2*dkdx3 + dkdx4) / 6;
        B(:, :, i) = (dkdu1 + 2*dkdu2 + 2*dkdu3 + dkdu4) / 6;
        d(:, i) = f - A(:, :, i)*x_i - B(:, :, i)*u_i;

    end
    
end

