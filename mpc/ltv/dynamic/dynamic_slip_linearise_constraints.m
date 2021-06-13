function [A, lb, ub] = dynamic_slip_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa)
%LINEARISE_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    ac_max = 6.5330;
    al_max = 10.0; 
    lr = 0.6183;    
    lf = 0.8672;

    [N_temp, N_x] = size(A_bar);
    [N_u, ~] = size(u_lin);
    N_steps = N_temp / N_x;

    % Friction ellipse constraint
    % g(x) = v^2 * delta / (lr + lf)
    % -5.0 <= g(x) <= 5.0
    
    C_bar = zeros(N_steps*2, N_steps*N_x);
    D_bar = zeros(N_steps*2, N_steps*N_u);
    g_bar = zeros(N_steps*2, 1);
    
    for i = 1:N_steps
        x = x_lin(:, i);
        u = u_lin(:, i);
        
        [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d, vf, denom_vf2] = A_curv_dyn(x, u, kappa);
        
        g0 = [-atan(vr);
              x(7) - atan(vf)];
        
        C = [0, 0, 0, denom_vr2*vr*x_d_hat_d/x_d_hat, -denom_vr2/x_d_hat, denom_vr2*lr/x_d_hat, 0;
             
             0, 0, 0, denom_vf2*vf*x_d_hat_d/x_d_hat, -denom_vf2/x_d_hat, -denom_vf2*lf/x_d_hat, 1;];
         
        D = [0, 0;
            0, 0];

        C_bar(2*i-1:i*2, (i-1)*N_x + 1 : i*N_x) = C;
        D_bar(2*i-1:i*2, (i-1)*N_u + 1 : i*N_u) = D;
        g_bar(2*i-1:i*2) = g0;
    end
        
    A = C_bar * B_bar + [D_bar, zeros(N_steps*2, 4)];
    
    const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:)) - D_bar*u_lin(:);
    lb = repmat([-0.1; -0.1], N_steps, 1) - const;
    ub = repmat([0.1; 0.1], N_steps, 1) - const;
    

end

% function [A, lb, ub] = dynamic_tyre_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa)
% %LINEARISE_CONSTRAINTS Summary of this function goes here
% %   Detailed explanation goes here
% 
%     ac_max = 6.5330;
%     al_max = 10.0; 
%     lr = 0.6183;    
%     lf = 0.8672;
% 
%     [N_temp, N_x] = size(A_bar);
%     [N_u, ~] = size(u_lin);
%     N_steps = N_temp / N_x;
% 
%     % Friction ellipse constraint
%     % g(x) = v^2 * delta / (lr + lf)
%     % -5.0 <= g(x) <= 5.0
%     
%     C_bar = zeros(N_steps*3, N_steps*N_x);
%     D_bar = zeros(N_steps*3, N_steps*N_u);
%     g_bar = zeros(N_steps*3, 1);
%     
%     for i = 1:N_steps
%         x = x_lin(:, i);
%         u = u_lin(:, i);
%         
%         [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d, vf, denom_vf2] = A_curv_dyn_lin(x, u, kappa);
%         
%         g0 = [(Fcr / (200*ac_max))^2 + (u(1) / al_max)^2;
%               -vr;
%               x(7) - vf];
%         
%         C = [0, 0, 0, ...
%             2*Fcr*Fcr_d*denom_vr2*vr*x_d_hat_d/x_d_hat / (200*ac_max)^2, ...
%              -2*Fcr*Fcr_d*denom_vr2/x_d_hat / (200*ac_max)^2, ...
%              2*Fcr*Fcr_d*denom_vr2*lr/x_d_hat / (200*ac_max)^2, ...
%              0;
%              
%              0, 0, 0, denom_vr2*vr*x_d_hat_d/x_d_hat, -denom_vr2/x_d_hat, denom_vr2*lr/x_d_hat, 0;
%              
%              0, 0, 0, denom_vf2*vf*x_d_hat_d/x_d_hat, -denom_vf2/x_d_hat, -denom_vf2*lf/x_d_hat, 1;];
%          
%         D = [2*u(1) / al_max^2, 0;
%             0, 0;
%             0, 0];
% 
%         C_bar(3*i-2:i*3, (i-1)*N_x + 1 : i*N_x) = C;
%         D_bar(3*i-2:i*3, (i-1)*N_u + 1 : i*N_u) = D;
%         g_bar(3*i-2:i*3) = g0;
%     end
%         
%     A = C_bar * B_bar + [D_bar, zeros(N_steps*3, 1)];
%     
%     const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:)) - D_bar*u_lin(:);
%     lb = repmat([-inf; -0.05; -0.05], N_steps, 1) - const;
%     ub = repmat([1.0; 0.05; 0.05], N_steps, 1) - const;
%     
% 
% end
% 


