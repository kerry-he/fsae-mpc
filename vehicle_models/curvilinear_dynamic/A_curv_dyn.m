function [A, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d, vf, denom_vf2] = A_curv_dyn(x, ~, kappa)
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
    m = 200;
    I = 200;
    lr = 0.6183;
    lf = 0.8672;
    
    g = 9.81;

    % Define states and controls
    s       = x(1);
    n       = x(2);
    mu      = x(3);
    x_d     = x(4);
    y_d     = x(5);
    theta_d = x(6);
    delta   = x(7); 
    
    x_d_hat = x_d + 5*exp(-x_d/5);
    x_d_hat_d = 1 - exp(-x_d/5);
    
    % Slip angles
    alpha_f = delta - atan((y_d + lf*theta_d) / x_d_hat);
    alpha_r = -atan((y_d - lr*theta_d) / x_d_hat);
    
    % Mass distribution
    Fzf = m*g * lf / (lr+lf);
    Fzr = m*g * lr / (lr+lf);
    
    % Pacejka magic formula
    B = 12.56;
    C = 1.38;
    D = 1.60;
    E = -0.58;
    
    Fcf = Fzf * D * sin(C * atan(B*alpha_f - E*(B*alpha_f - atan(B*alpha_f))));
    Fcr = Fzr * D * sin(C * atan(B*alpha_r - E*(B*alpha_r - atan(B*alpha_r))));
    
    Fcf_d = Fzf * D * cos(C * atan(B*alpha_f - E*(B*alpha_f - atan(B*alpha_f)))) ...
                * C / (1 + (B*alpha_f - E*(B*alpha_f - atan(B*alpha_f)))^2) ...
                * (B - E * (B - B / (1 + B^2 * alpha_f^2)));
            
    Fcr_d = Fzr * D * cos(C * atan(B*alpha_r - E*(B*alpha_r - atan(B*alpha_r)))) ...
                * C / (1 + (B*alpha_r - E*(B*alpha_r - atan(B*alpha_r)))^2) ...
                * (B - E * (B - B / (1 + B^2 * alpha_r^2)));            
    
	% Define common constants
    k = kappa(s);
    denom_nk = 1 / (1 - n * k); 
    vf = (y_d + lf*theta_d) / x_d_hat;
    vr = (y_d - lr*theta_d) / x_d_hat;
    denom_vf2 = 1 / (1 + vf^2);
    denom_vr2 = 1 / (1 + vr^2);
    
    % Partial derivatives
    s_n = (x_d * cos(mu) - y_d * sin(mu))*denom_nk^2 * k;
    s_mu = (-x_d * sin(mu) - y_d * cos(mu))*denom_nk;
    s_xd = cos(mu) * denom_nk;
    s_yd = -sin(mu) * denom_nk;
    
    n_mu = x_d * cos(mu) - y_d * sin(mu);
    n_xd = sin(mu);
    n_yd = cos(mu);
    
    mu_n = -s_n * k;
    mu_mu = -s_mu * k;
    mu_xd = -s_xd * k;
    mu_yd = -s_yd * k;
    mu_thetad = 1;
    
    xd_xd = -Fcf_d * denom_vf2 * vf * sin(delta) * x_d_hat_d / (m * x_d_hat);
    xd_yd = (Fcf_d * denom_vf2 * sin(delta) / x_d_hat + m * theta_d) / m;
    xd_thetad = (Fcf_d * denom_vf2 * lf  * sin(delta) / x_d_hat + m * y_d) / m;
    xd_delta = (-Fcf * cos(delta) - Fcf_d * sin(delta)) / m;
    
    yd_xd = (Fcr_d * denom_vr2 * vr * x_d_hat_d / x_d_hat + Fcf_d * denom_vf2 * vf * cos(delta) * x_d_hat_d / x_d_hat - m * theta_d) / m;
    yd_yd = (-Fcr_d  * denom_vr2 / x_d_hat - Fcf_d * denom_vf2 / x_d_hat * cos(delta)) / m;
    yd_thetad = (Fcr_d  * denom_vr2 * lr / x_d_hat - Fcf_d * denom_vf2 * lf / x_d_hat * cos(delta) - m * x_d_hat) / m;
    yd_delta = (-Fcf * sin(delta) + Fcf_d * cos(delta)) / m;
    
    t_xd = (lf * Fcf_d * denom_vf2 * vf * cos(delta) * x_d_hat_d / x_d_hat - lr * Fcr_d * denom_vr2 * vr * x_d_hat_d / x_d_hat) / I;
    t_yd = (-lf * Fcf_d * denom_vf2 * cos(delta) / x_d_hat + lr * Fcr_d * denom_vr2 / x_d_hat) / I;
    t_thetad = (-lf * Fcf_d * denom_vf2 * lf * cos(delta) / x_d_hat - lr * Fcr_d * denom_vr2 * lr / x_d_hat) / I;
    t_delta = (-lf * Fcf * sin(delta) + lf * Fcf_d * cos(delta)) / I;
    
    % Populate matrices
    A = [0          s_n        s_mu        s_xd       s_yd       0          0;
         0          0          n_mu        n_xd       n_yd       0          0;
         0          mu_n       mu_mu       mu_xd      mu_yd      mu_thetad  0;
         0          0          0           xd_xd      xd_yd      xd_thetad  xd_delta;
         0          0          0           yd_xd      yd_yd      yd_thetad  yd_delta;
         0          0          0           t_xd       t_yd       t_thetad   t_delta;
         0          0          0           0          0          0          0];

end

