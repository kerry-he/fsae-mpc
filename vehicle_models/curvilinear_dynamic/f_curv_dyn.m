function f = f_curv_dyn(x, u, kappa)
%F_CART_DYN The dynamic equations of the kinematic bicycle model in
%Cartesian coodinates.
%   INPUT:
%       x - State vector to linearise at [s; n; mu; x_d, y_d, theta_d, delta]
%       u - Control vector to linearise at[Fx; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);  
%   OUTPUT:
%       f - Dynamic model equations

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
    x_d     = max(x(4), 0.01);
    y_d     = x(5);
    theta_d = x(6);
    delta   = x(7);
    
    Fx      = u(1)*m;
    delta_d = u(2);       
    
	% Define common constants
    k = kappa(s);
    denom_nk = 1 / (1 - n * k); 
    
    % Slip angles
    alpha_f = delta - atan((y_d + lf*theta_d) / x_d);
    alpha_r = -atan((y_d - lr*theta_d) / x_d);
    
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
    
    % Populate matrix
    f = [(x_d * cos(mu) - y_d * sin(mu))*denom_nk;
         x_d * sin(mu) + y_d * cos(mu);
         theta_d - (x_d * cos(mu) - y_d * sin(mu))*denom_nk*k;
         (Fx - Fcf*sin(delta) + m*y_d*theta_d) / m;
         (Fcr + Fcf*cos(delta) - m*x_d*theta_d) / m;
         (lf*Fcf*cos(delta) - lr*Fcr) / I;
         delta_d];

end

