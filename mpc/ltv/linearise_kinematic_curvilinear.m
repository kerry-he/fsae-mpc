function [A, B, d] = linearise_kinematic_curvilinear(x, u, kappa)
%LINEARISE_KINEMATIC_CURVILINEAR Linearises the dynamics of the kinematic
%bicycle model using a curvilinear coordinate frame at a given setpoint
%   INPUTS:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);   
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
    
    for i = 1:N_steps
        % Precompute common expressions
        k = kappa(x(1, i));
        beta = atan(lr_ratio * tan(x(5, i)));
        s_mu_beta = sin(x(3, i) + beta);
        c_mu_beta = cos(x(3, i) + beta);
        beta_d = lr_ratio * sec(x(5, i))^2 / ...
            (1 + (lr_ratio * tan(x(5, i)))^2);
        denom_nk = 1 / (1 - x(2, i) * k); 

        % Compute partial derivatives
        s_n = x(4, i)*c_mu_beta * denom_nk^2 * k;
        s_mu = -x(4, i)*s_mu_beta * denom_nk;
        s_v = c_mu_beta * denom_nk;
        s_delta = -x(4, i)*s_mu_beta * denom_nk * beta_d;

        n_mu = x(4, i)*c_mu_beta;
        n_v = s_mu_beta;
        n_delta = x(4, i)*c_mu_beta * beta_d;

        mu_n = -s_n * k;
        mu_mu = -s_mu * k;
        mu_v = sin(beta)/lr - s_v * k;
        mu_delta = x(4, i)*cos(beta)*beta_d/lr - s_delta * k;

        % Populate matrices
        A(:, :, i) = sparse([0   s_n    s_mu    s_v    s_delta;
                             0   0      n_mu    n_v    n_delta;
                             0   mu_n   mu_mu   mu_v   mu_delta;
                             0   0      0       0      0;
                             0   0      0       0      0]);

        B(:, :, i) = sparse([0 0;
                             0 0;
                             0 0;
                             1 0;
                             0 1]);

        d(:, i) = [x(4, i)*c_mu_beta*denom_nk;
                   x(4, i)*s_mu_beta;
                   x(4, i)*sin(beta)/lr - x(4, i)*c_mu_beta*denom_nk*k;
                   u(1, i);
                   u(2, i)] - A(:, :, i)*x(:, i) - B(:, :, i)*u(:, i);
    end
    
end

