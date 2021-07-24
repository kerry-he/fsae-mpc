function [A, lb, ub] = dynamic_tyre_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa)
%LINEARISE_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    ac_max = 9.163;
    al_max = 10.0; 
    lr = 0.6183;    

    [N_temp, N_x] = size(A_bar);
    [N_u, ~] = size(u_lin);
    N_steps = N_temp / N_x;

    % Friction ellipse constraint
    % g(x) = (ac / ac_max)^2 + (al / al_max)^2
    % g(x) <= 1.0
    
    % Linearise tyre constraints
    N = 12;
    theta = linspace(0, 2*pi, N + 1);
    ac_list = ac_max * sin(theta);
    al_list = al_max * cos(theta);
    dac = ac_list(2:end)-ac_list(1:N);
    dal = al_list(2:end)-al_list(1:N);
    
    
    C_bar = zeros(N_steps*N, N_steps*N_x);
    D_bar = zeros(N_steps*N, N_steps*N_u);
    g_bar = zeros(N_steps*N, 1);
    
    for i = 1:N_steps
        x = x_lin(:, i);
        u = u_lin(:, i);
        
        [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d] = A_curv_dyn(x, u, kappa);
        
        
        g0 = zeros(N, 1);
        C = zeros(N, N_x);
        D = zeros(N, N_u);
        for j = 1:N
            g0(j) = (u(1) - al_list(j))*dac(j) - (Fcr/280 - ac_list(j))*dal(j);
            
            C(j, :) = [0, 0, 0, ...
                       -dal(j)*Fcr_d*denom_vr2*vr*x_d_hat_d/x_d_hat / 280, ...
                       dal(j)*Fcr_d*denom_vr2/x_d_hat / 280, ...
                       -dal(j)*Fcr_d*denom_vr2*lr/x_d_hat / 280, ...
                       0];
                   
            D(j, :) = [dac(j), 0];                   
        end

        C_bar(N*(i-1)+1:i*N, (i-1)*N_x + 1 : i*N_x) = C;
        D_bar(N*(i-1)+1:i*N, (i-1)*N_u + 1 : i*N_u) = D;
        g_bar(N*(i-1)+1:i*N) = g0;
    end
        
    A = C_bar * B_bar + [D_bar, zeros(N_steps*N, 4)];
    
    const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:)) - D_bar*u_lin(:);
    lb = -inf(N_steps*N, 1);
    ub = zeros(N_steps*N, 1) - const;
    

end