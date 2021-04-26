function f = f_cart_kin(x, u)
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
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);
    
    beta = atan(lr_ratio * tan(x(5)));
    
    % Populate matrix
    f = [x(4) * cos(x(3) + beta);
         x(4) * sin(x(3) + beta);
         x(4) / lr * sin(beta);
         u(1);
         u(2)];

end

