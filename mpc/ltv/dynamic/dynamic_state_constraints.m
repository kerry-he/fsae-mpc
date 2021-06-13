function [B_bar, xA, lbA, ubA] = dynamic_state_constraints(A_bar, B_bar, d_bar, x0, lb, ub, state_idx, soft_idx, x_lin, u_lin, kappa)
%SOFT_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    [N_temp, N_x] = size(A_bar);
    N_steps = N_temp / N_x;
    N_state = length(state_idx);
    N_soft = length(soft_idx);

    % Append B_bar matrix to include soft constraints
    B_bar = [B_bar, zeros(N_x*N_steps, 4)];
    
    % Define state constraints
    state_idx_full = zeros(N_state*N_steps, 1);
    soft_idx_full = zeros(N_soft*N_steps, 1);
    for i = 1:N_state
        state_idx_full((i - 1)*N_steps+1:i*N_steps) = state_idx(i):N_x:N_x*N_steps;
    end
    for i = 1:N_soft
        soft_idx_full((i - 1)*N_steps+1:i*N_steps) = soft_idx(i):N_x:N_x*N_steps;
    end
    constraint_idx = [state_idx_full; soft_idx_full];
    
    % Segment out parts of matrix which correspond to variables which
    % constraints are applied to
    A = A_bar(constraint_idx, :);
    xA = B_bar([constraint_idx; soft_idx_full], :);
    D = d_bar(constraint_idx, :);
    const = A*x0 + D;
    
    
    
    % Define lower and upper bounds
    
    lbA = lb - const;
    ubA = ub - const;
    
    lbA = [lbA; -ones(N_soft*N_steps, 1)*1e10];
    ubA = [ubA(1:N_state*N_steps); ones(N_soft*N_steps, 1)*1e10; ubA(N_state*N_steps+1:end)];
    
    % Modify B_bar matrix to account for soft constraints
    xA(end-2*N_soft*N_steps+1:end, end-3) = [ones(N_soft*N_steps, 1); -ones(N_soft*N_steps, 1)];
    
    
    % Functional constraints
    [A_ay, lb_ay, ub_ay] = dynamic_slip_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa);
    xA = [xA; A_ay; A_ay];
    lbA = [lbA; lb_ay; -inf*ones(N_steps*2, 1)];
    ubA = [ubA; inf*ones(N_steps*2, 1); ub_ay];
    xA(end-4*N_steps+1:end, end-2:end-1) = [repmat([1, 0; 0, 1], N_steps, 1); -repmat([1, 0; 0, 1], N_steps, 1)];   
    
    
    [A_ay, lb_ay, ub_ay] = dynamic_tyre_linearise_constraints(A_bar, B_bar, d_bar, x_lin, u_lin, x0, kappa);
    xA = [xA; A_ay];
    lbA = [lbA; lb_ay];
    ubA = [ubA; ub_ay];
    xA(end-N_steps+1:end, end) = -ones(N_steps, 1);    
end

