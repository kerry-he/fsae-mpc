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

    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);

    % Preallocate
    A = zeros(N_x, N_x, N_steps);
    B = zeros(N_x, N_u, N_steps);
    d = zeros(N_x, N_steps);
    
    rk_factor = [1 2 2 1];
    for i = 1:N_steps
        x_k = x(:, i);
        for j = 1:4
            % Precompute common expressions
            k = kappa(x_k(1));
            beta = atan(lr_ratio * tan(x_k(5)));
            s_mu_beta = sin(x_k(3) + beta);
            c_mu_beta = cos(x_k(3) + beta);
            beta_d = lr_ratio * sec(x_k(5))^2 / ...
                (1 + (lr_ratio * tan(x_k(5)))^2);
            denom_nk = 1 / (1 - x_k(2) * k); 

            % Compute partial derivatives
            s_n = x_k(4)*c_mu_beta * denom_nk^2 * k;
            s_mu = -x_k(4)*s_mu_beta * denom_nk;
            s_v = c_mu_beta * denom_nk;
            s_delta = -x_k(4)*s_mu_beta * denom_nk * beta_d;

            n_mu = x_k(4)*c_mu_beta;
            n_v = s_mu_beta;
            n_delta = x_k(4)*c_mu_beta * beta_d;

            mu_n = -s_n * k;
            mu_mu = -s_mu * k;
            mu_v = sin(beta)/lr - s_v * k;
            mu_delta = x_k(4)*cos(beta)*beta_d/lr - s_delta * k;

            % Populate matrices
            A_k = [0      s_n    s_mu    s_v    s_delta;
                   0      0      n_mu    n_v    n_delta;
                   0      mu_n   mu_mu   mu_v   mu_delta;
                   0      0      0       0      0;
                   0      0      0       0      0];
            A(:, :, i) = A(:, :, i) + A_k*rk_factor(j);

            B_k = [0 0;
                   0 0;
                   0 0;
                   1 0;
                   0 1];
            B(:, :, i) = B(:, :, i) + B_k*rk_factor(j);

            f = [x_k(4)*c_mu_beta*denom_nk;
                 x_k(4)*s_mu_beta;
                 x_k(4)*sin(beta)/lr - x_k(4)*c_mu_beta*denom_nk*k;
                 u(1, i);
                 u(2, i)];

            d(:, i) = d(:, i) + (f - A_k*x_k - B_k*u(:, i)) * rk_factor(j);

            x_k = x(:, i) + f*dt / rk_factor(mod(j, 4)+1);
        end
    end
    
    A = A / 6;
    B = B / 6;
    d = d / 6;
    
end

