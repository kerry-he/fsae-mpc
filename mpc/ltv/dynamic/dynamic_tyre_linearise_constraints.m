function [A, lb, ub] = dynamic_tyre_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa)
%LINEARISE_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    ac_max = 6.5330;
    al_max = 10.0; 
    lr = 0.6183;    

    [N_temp, N_x] = size(A_bar);
    [N_u, ~] = size(u_lin);
    N_steps = N_temp / N_x;

    % Friction ellipse constraint
    % g(x) = v^2 * delta / (lr + lf)
    % -5.0 <= g(x) <= 5.0
    
    C_bar = zeros(N_steps, N_steps*N_x);
    D_bar = zeros(N_steps, N_steps*N_u);
    g_bar = zeros(N_steps, 1);
    
    for i = 1:N_steps
        x = x_lin(:, i);
        u = u_lin(:, i);
        
        [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d] = A_curv_dyn_lin(x, u, kappa);
        
        g0 = (Fcr / (200*ac_max))^2 + (u(1) / al_max)^2;
        
        C = [0, 0, 0, ...
            2*Fcr*Fcr_d*denom_vr2*vr*x_d_hat_d/x_d_hat / (200*ac_max)^2, ...
             -2*Fcr*Fcr_d*denom_vr2/x_d_hat / (200*ac_max)^2, ...
             2*Fcr*Fcr_d*denom_vr2*lr/x_d_hat / (200*ac_max)^2, ...
             0;];
         
        D = [2*u(1) / al_max^2, 0];

        C_bar(i, (i-1)*N_x + 1 : i*N_x) = C;
        D_bar(i, (i-1)*N_u + 1 : i*N_u) = D;
        g_bar(i) = g0;
    end
        
    A = C_bar * B_bar + [D_bar, zeros(N_steps, 4)];
    
    const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:)) - D_bar*u_lin(:);
    lb = -inf(N_steps, 1) - const;
    ub = ones(N_steps, 1) - const;
    

end





