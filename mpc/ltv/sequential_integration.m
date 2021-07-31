function [A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt)
%MAKE_SEQUENTIAL Sequential integration of a given linearised system
%   INPUTS:
%       Parameters of the linearised differential equation: 
%       dx/dt = Ax + Bu + d
%       dt - Integration time step
%   OUTPUTS:
%       Parameters of the sequentially integrated system:
%           x = A_bar*x0 + B_bar*u + d_bar

    % Define constant sizes
    [N_x, N_steps] = size(d);
    [~, N_u, ~] = size(B); 

    % Euler integration
    A = A*dt + repmat(eye(N_x), 1, 1, N_steps);
    B = B*dt;
    d = d*dt;
    
    % Construct full matrices
    A_bar = zeros(N_x*N_steps, N_x);
    A_bar(1:N_x, :) = A(:, :, 1);
    for i = 2:N_steps
        A_bar((i - 1)*N_x + 1:i*N_x, :) = ...
        	A(:, :, i) * A_bar((i - 2)*N_x + 1:(i - 1)*N_x, :);
    end
    
    B_bar = zeros(N_x*N_steps, N_u*N_steps);
    for i = 1:N_steps
        B_bar((i - 1)*N_x + 1:i*N_x, (i - 1)*N_u + 1:i*N_u) = B(:, :, i);
        
        for j = i + 1:N_steps
            B_bar((j - 1)*N_x + 1:j*N_x, (i - 1)*N_u + 1:i*N_u) = A(:, :, j) * ...
                B_bar((j - 2)*N_x + 1:(j - 1)*N_x, (i - 1)*N_u + 1:i*N_u);
        end
    end
    
    D = zeros(N_x*N_steps);
    for i = 1:N_steps
        D((i - 1)*N_x + 1:i*N_x, (i - 1)*N_x + 1:i*N_x) = eye(N_x);
        
        for j = i + 1:N_steps
            D((j - 1)*N_x + 1:j*N_x, (i - 1)*N_x + 1:i*N_x) = A(:, :, j) * ...
                D((j - 2)*N_x + 1:(j - 1)*N_x, (i - 1)*N_x + 1:i*N_x);
        end
    end    
    d_bar = D*d(:);

end

