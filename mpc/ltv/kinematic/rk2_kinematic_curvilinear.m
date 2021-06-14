function [A, B, d] = rk2_kinematic_curvilinear(x, u, kappa, dt)
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
        k1 = f_curv_kin(x_i, u_i, kappa);
        k2 = f_curv_kin(x_i + k1*dt / 2, u_i, kappa);
        f = k2;
                
        % Calculate RK4 partial derivatives w.r.t. states 
        dfdx1 = A_curv_kin(x_i, u_i, kappa);
        dfdx2 = A_curv_kin(x_i + k1*dt / 2, u_i, kappa);
        
        dkdx1 = dfdx1;
        dkdx2 = dfdx2 * (I + dkdx1*dt / 2);
        
        % Calculate RK4 partial derivatives w.r.t. controls
        dkdu1 = B_curv_kin(x_i, u_i, kappa);
        dkdu2 = B_curv_kin(x_i + k1*dt / 2, u_i, kappa) + dfdx2*dkdu1*dt / 2;
        
        % Calculate final matrices
        A(:, :, i) = dkdx2;
        B(:, :, i) = dkdu2;
        d(:, i) = f - A(:, :, i)*x_i - B(:, :, i)*u_i;

    end
    
end

