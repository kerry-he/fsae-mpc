function [A, lb, ub] = linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa)
%LINEARISE_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    ac_max = 9.1630;
    al_max = 10.0; 
    lr = 0.6183;    
    lf = 0.8672;
    K_vel = 1.6;
    K_steer = 5.0;

    [N_temp, N_x] = size(A_bar);
    [N_u, ~] = size(u_lin);
    N_steps = N_temp / N_x;

    % Friction ellipse constraint
    % g(x) = (ac / ac_max)^2 + (al / al_max)^2
    % g(x) <= 1.0
    
    N_other = 7;
    
    % Linearise tyre constraints
    N_tyre = 8;
    theta = linspace(0, 2*pi, N_tyre + 1);
    ac_list = ac_max * sin(theta);
    al_list = al_max * cos(theta);
    dac = ac_list(2:end)-ac_list(1:N_tyre);
    dal = al_list(2:end)-al_list(1:N_tyre);
    
    
    N = N_other + N_tyre;
    
    C_bar = zeros(N_steps*N, N_steps*N_x);
    D_bar = zeros(N_steps*N, N_steps*N_u);
    g_bar = zeros(N_steps*N, 1);
    
    
    
    for i = 1:N_steps
        x = x_lin(:, i);
        u = u_lin(:, i);
        
        [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d, vf, denom_vf2] = A_curv_dyn(x, u, kappa);

        g0 = [K_steer * (u(2) - x(7));
              x(2);
              x(2);
              -atan(vr);
              x(7) - atan(vf);
              -atan(vr);
              x(7) - atan(vf)];            

        C = [0, 0, 0, 0, 0, 0, -K_steer;
             0, 1, 0, 0, 0, 0, 0;
             0, 1, 0, 0, 0, 0, 0;
             0, 0, 0, denom_vr2*vr*x_d_hat_d/x_d_hat, -denom_vr2/x_d_hat, denom_vr2*lr/x_d_hat, 0;
             0, 0, 0, denom_vf2*vf*x_d_hat_d/x_d_hat, -denom_vf2/x_d_hat, -denom_vf2*lf/x_d_hat, 1;
             0, 0, 0, denom_vr2*vr*x_d_hat_d/x_d_hat, -denom_vr2/x_d_hat, denom_vr2*lr/x_d_hat, 0;
             0, 0, 0, denom_vf2*vf*x_d_hat_d/x_d_hat, -denom_vf2/x_d_hat, -denom_vf2*lf/x_d_hat, 1];

        D = [0, K_steer;
             0, 0;
             0, 0;
             0, 0;
             0, 0;
             0, 0;
             0, 0];
         
        for j = 1:N_tyre
            g0(N_other + j) = (K_vel*(u(1) - x(4)) - al_list(j))*dac(j) - (Fcr/280 - ac_list(j))*dal(j);
            
            C(N_other + j, :) = [0, 0, 0, ...
                                 -dal(j)*Fcr_d*denom_vr2*vr*x_d_hat_d/x_d_hat / 280 - K_vel*dac(j), ...
                                 dal(j)*Fcr_d*denom_vr2/x_d_hat / 280, ...
                                 -dal(j)*Fcr_d*denom_vr2*lr/x_d_hat / 280, ...
                                 0];
                   
            D(N_other + j, :) = [K_vel*dac(j), 0];
        end    

        C_bar(N*(i-1)+1:i*N, (i-1)*N_x + 1 : i*N_x) = C;
        D_bar(N*(i-1)+1:i*N, (i-1)*N_u + 1 : i*N_u) = D;
        g_bar(N*(i-1)+1:i*N) = g0;
    end
        
    A = C_bar * B_bar + [D_bar, repmat([0,  0,  0,  0;
                                       -1,  0,  0,  0;
                                        1,  0,  0,  0;
                                        0, -1,  0,  0;
                                        0,  0, -1,  0;
                                        0,  1,  0,  0;
                                        0,  0,  1,  0;
                                        repmat([0, 0, 0, 1], N_tyre, 1)], N_steps, 1);];
    
    const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:)) - D_bar*u_lin(:);
    lb = repmat([-0.8; -inf;  -0.75; -inf; -inf; -0.1; -0.1; -inf(N_tyre, 1)], N_steps, 1) - const;
    ub = repmat([0.8;   0.75;  inf;   0.1;  0.1;  inf;  inf; ones(N_tyre, 1)], N_steps, 1) - const;
    

end