function [H, f] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R, R_soft)
%GENERATE_QP Creates quadratic programming parameters for a given LTV
%system
%   INPUTS:
%       Parameters of the sequentially integrated system:
%           x = A_bar*x0 + B_bar*u + d_bar
%       x_0 - Initial conditions
%       x_ref - State reference vector
%       Q - State weight vector
%       Q_terminal - Terminal state weight vector
%       R - Control weight vector
%       R_soft - Soft constraints weight vector
%   OUTPUTS: 
%       Parameters of the QP: 
%           1/2x'Hx + f'x

    % Define constant sizes
    [N_x, N_steps] = size(x_ref);
    N_u = length(R);
    N_soft = length(R_soft);

    % Construct full weight matrices
    Q_bar = spdiags([repmat(Q(:), N_steps - 1, 1); Q_terminal(:)], 0, ...
        N_steps*N_x, N_steps*N_x);
    R_bar = spdiags([repmat(R(:), N_steps, 1); repmat(R_soft(:), N_steps, 1)], ...
        0, N_steps*(N_u+N_soft), N_steps*(N_u+N_soft));
    
    % Define QP parameters
    H = 2 * (B_bar' * Q_bar * B_bar + R_bar);
    f = 2 * B_bar' * Q_bar * (A_bar * x0 + d_bar - x_ref(:));
    
end

