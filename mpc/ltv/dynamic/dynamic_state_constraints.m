function [B_bar, xA, lbA, ubA] = dynamic_state_constraints(A_bar, B_bar, d_bar, x0, lb, ub, state_idx, soft_idx, x_lin, u_lin, kappa)
%SOFT_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    [N_temp, N_x] = size(A_bar);
    N_steps = N_temp / N_x;
    N_state = length(state_idx);
    N_soft = length(soft_idx);

    % Append B_bar matrix to include soft constraints
    B_bar = [B_bar, zeros(N_x*N_steps, 4)];
    
    % Functional constraints
    [xA, lbA, ubA] = linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa);
end

