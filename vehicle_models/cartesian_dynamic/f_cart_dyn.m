function f = f_cart_dyn(x, u)
%F_CART_KIN The dynamic equations of the kinematic bicycle model in
%Cartesian coodinates.
%   INPUT:
%       x - State vector to linearise at [s; n; mu; v; delta]
%       u - Control vector to linearise at[a; delta_d]
%       kappa - Anonymous function of spline curvature:
%       	kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl);  
%   OUTPUT:
%       f - Dynamic model equations

    % Define vehicle constants
    m = 280;
    I = 230;
    lr = 0.6183;
    lf = 0.8672;
    
    g = 9.81;
    
    % Define states and controls
    theta   = x(3);
    x_d     = x(4);
    y_d     = x(5);
    theta_d = x(6);
    delta   = x(7);
    
    Fx      = u(1);
    delta_d = u(2);    
    
    % Slip angles
    alpha_f = delta - atan((y_d + lf*theta_d) / (x_d + 0.01));
    alpha_r = -atan((y_d - lr*theta_d) / (x_d + 0.01));
    
    % Mass distribution
    Fzf = m*g * lr / (lr+lf);
    Fzr = m*g * lf / (lr+lf);
    
    % Pacejka magic formula
    B = 12.56;
    C = 1.38;
    D = 1.60;
    E = -0.58;
    
    Fcf = Fzf * D * sin(C * atan(B*alpha_f - E*(B*alpha_f - atan(B*alpha_f))));
    Fcr = Fzr * D * sin(C * atan(B*alpha_r - E*(B*alpha_r - atan(B*alpha_r))));
    
    % Populate matrix
    f = [x_d * cos(theta) - y_d * sin(theta);
         x_d * sin(theta) + y_d * cos(theta);
         theta_d;
         (Fx - Fcf*sin(delta) + m*y_d*theta_d) / m;
         (Fcr + Fcf*cos(delta) - m*x_d*theta_d) / m;
         (lf*Fcf*cos(delta) - lr*Fcr) / I;
         delta_d];

end

