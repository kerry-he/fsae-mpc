function [A, lb, ub] = linearise_constraints(A_bar, B_bar, d_bar, x_lin, x0)
%LINEARISE_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    lr = 0.6183;
    lf = 0.8672;

    [N_temp, N_x] = size(A_bar);
    N_steps = N_temp / N_x;

    % Friction ellipse constraint
    % g(x) = v^2 * delta / (lr + lf)
    % -5.0 <= g(x) <= 5.0
    
    C_bar = zeros(N_steps, N_steps*N_x);
    g_bar = zeros(N_steps, 1);
    
    for i = 1:N_steps
        x = x_lin(:, i);
        
        g0 = x(4)^2 * x(5) / (lr + lf);
        C = [0, 0, 0, 2*x(4)*x(5), x(4)^2] / (lf + lr);
        
        C_bar(i, (i-1)*N_x + 1 : i*N_x) = C;
        g_bar(i) = g0;
    end
        
    A = C_bar * B_bar;
    
    const = g_bar + C_bar * (A_bar*x0 + d_bar - x_lin(:));
    lb = -5.0 - const;
    ub = 5.0 - const;
    

end

