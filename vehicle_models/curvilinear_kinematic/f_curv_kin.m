function f = f_curv_kin(x, u, kappa)
%F_CURV_KIN The dynamic equations of the kinematic bicycle model in
%curvilinear coodinates.
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

	% Define common constants
    k = kappa(x(1));
    beta = atan(lr_ratio * tan(x(5)));
    s_mu_beta = sin(x(3) + beta);
    c_mu_beta = cos(x(3) + beta);
    denom_nk = 1 / (1 - x(2) * k); 

    % Populate matrix
    f = [x(4)*c_mu_beta*denom_nk;
         x(4)*s_mu_beta;
         x(4)*sin(beta)/lr - x(4)*c_mu_beta*denom_nk*k;
         u(1);
         u(2)];

end

