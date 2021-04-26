function A = A_curv_kin(x, ~, kappa, kappa_d)
%F_CURV_KIN Calculates the Jacobian of f_curv_kin with respect to the state
%variables
%   INPUT:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);  
%       kappa_d - (Optional) Anonymous function of spline curvature
%           derivative
%   OUTPUT:
%       A - Jacobian of vehicle model with respect to state variables

    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);

	% Define common constants
    k = kappa(x(1));
    beta = atan(lr_ratio * tan(x(5)));
    s_mu_beta = sin(x(3) + beta);
    c_mu_beta = cos(x(3) + beta);
    beta_d = lr_ratio * sec(x(5))^2 / (1 + (lr_ratio * tan(x(5)))^2);
    denom_nk = 1 / (1 - x(2) * k); 

    % Compute partial derivatives
    s_s = 0;
    s_n = x(4)*c_mu_beta * denom_nk^2 * k;
    s_mu = -x(4)*s_mu_beta * denom_nk;
    s_v = c_mu_beta * denom_nk;
    s_delta = -x(4)*s_mu_beta * denom_nk * beta_d;

    n_mu = x(4)*c_mu_beta;
    n_v = s_mu_beta;
    n_delta = x(4)*c_mu_beta * beta_d;

    mu_s = 0;
    mu_n = -s_n * k;
    mu_mu = -s_mu * k;
    mu_v = sin(beta)/lr - s_v * k;
    mu_delta = x(4)*cos(beta)*beta_d/lr - s_delta * k;
    
    if exist('kappa_d','var')
        k_d = kappa_d(x(1));
        s_s = x(4)*c_mu_beta * denom_nk^2 * k_d * x(2);
        mu_s = -x(4)*c_mu_beta * denom_nk * k_d - s_s * k;
    end

    % Populate matrices
    A = [s_s    s_n    s_mu    s_v    s_delta;
         0      0      n_mu    n_v    n_delta;
         mu_s   mu_n   mu_mu   mu_v   mu_delta;
         0      0      0       0      0;
         0      0      0       0      0];

end

