function [B_bar, xA, lbA, ubA] = state_constraints(A_bar, B_bar, d_bar, x0, lb, ub, state_idx, soft_idx)
%SOFT_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

    [N_temp, N_x] = size(A_bar);
    N_steps = N_temp / N_x;
    N_state = length(state_idx);
    N_soft = length(soft_idx);

    % Append B_bar matrix to include soft constraints
    B_bar = [B_bar, zeros(N_x*N_steps, N_soft*N_steps)];
    
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
    xA = B_bar(constraint_idx, :);
    D = d_bar(constraint_idx, :);
    const = A*x0 + D;
    
    % Define lower and upper bounds
    lbA = lb - const;
    ubA = ub - const;
    
    % Modify B_bar matrix to account for soft constraints
    xA(end-N_soft*N_steps+1 : end, end-N_soft*N_steps+1 : end) = eye(N_soft*N_steps);    

end

