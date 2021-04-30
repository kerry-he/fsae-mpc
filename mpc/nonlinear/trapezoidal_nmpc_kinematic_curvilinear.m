function [x, info] = trapezoidal_nmpc_kinematic_curvilinear(x0, x_ref, kappa, kappa_d, dt, x_init, info)
%NMPC_KINMATIC_CURVILINEAR Computes a NMPC step for a kinematic bicycle
%model using a curvilinear coordinate frame.
%   INPUTS:
%       x0 - Initial state
%       x_ref - Reference trajectory [x_ref_1, x_ref_2, ...]
%       kappa - Spline function
%       dt - Time step
%       info - IPOPT information from previous iteration
%   OUTPUTS:
%       x - Optimised variable
%       info - Information about solved problem

    % Warm starting
    if length(info) == 1
        options.ipopt.warm_start_init_point = 'yes';
        options.zl = info.zl;
        options.zu = info.zu; 
        options.lambda = info.lambda;
    end
    

    % Define constants
    [N_x, N_steps] = size(x_ref);
    N_steps = N_steps + 1;
    N_u = 2;
    
    x_ref = [x0, x_ref];
    x_ref = [x_ref; zeros(N_u, N_steps)];
    

    % Defining cost weights
    Q = [5; 500; 2000; 0; 0];
    Q_terminal = Q * 10;
    R = [10, 10];
    
    Q_vec = [repmat([Q(:); R(:)], N_steps - 1, 1); Q_terminal(:); R(:)];
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { x0, x_ref, kappa, kappa_d, Q_bar, N_x, N_u, N_steps, dt };

    % The constraint functions are bounded from below by zero.
    options.lb = repmat([-inf; -inf; -inf; 0; -0.4; -10.0; -0.4], N_steps, 1); % Lower bound on optimization variable
    options.ub = repmat([inf; inf; inf; inf; 0.4; 10.0; 0.4], N_steps, 1); % Upper bound on optimization variable
    options.cl = zeros(N_x*N_steps, 1); % Lower bound on constraint function
    options.cu = zeros(N_x*N_steps, 1); % Upper bound on constraint function
    
    % Set IPOPT options
    options.ipopt.print_level           = 0;
    options.ipopt.max_iter              = 5000;
    options.ipopt.tol                   = 1e-6; % OR 1e-5
    options.ipopt.hessian_approximation = 'limited-memory';
%     options.ipopt.derivative_test       = 'first-order';
%     options.ipopt.derivative_test_tol   = 1e-2;
    
    % Define callback functions
    funcs.objective         = @objective; % Objective function (Required)
    funcs.gradient          = @gradient; % Gradient of objective (Required)
    
    funcs.constraints       = @constraints; % Constraint function, corresponds to cu/cl (Optional)
    funcs.jacobian          = @jacobian; % Jacobian of the constraints at the current points (Optional)
    funcs.jacobianstructure = @jacobianstructure; %Structure of Jacobian (Optional)

    % Run IPOPT.
    x_init(1:end-(N_x+N_u)) = x_init(N_x+N_u+1:end);
    x_init(end-(N_x+N_u)+1 : end-(N_u)) = x_init(end-(N_x+N_u)+1 : end-(N_u))...
        + dt*f_curv_kin(x_init(end-(N_x+N_u)+1 : end-(N_u)), x_init(end-N_u+1 : end), kappa);    
    [x, info] = ipopt_auxdata(x_init(:), funcs, options);  
    
% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [~, x_ref, ~, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x - x_ref(:);
    f = x_error' * Q_bar * x_error;

% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [~, x_ref, ~, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x - x_ref(:);
    g = 2 * Q_bar * x_error;

% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [x0, ~, kappa, ~, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    % Preallocate
    c = zeros(N_x*N_steps, 1);
    
    c(1:N_x) = x(1:N_x) - x0;
    
    for i = 1:N_steps-1
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        x_i_1 = x(i*(N_x+N_u) + 1:i*(N_x+N_u) + N_x);
        u_i_1 = x(i*(N_x+N_u) + N_x + 1  :  (i+1)*(N_x+N_u));
        
        % Calculate dynamic model
        f_i = f_curv_kin(x_i, u_i, kappa);
        f_i_1 = f_curv_kin(x_i_1, u_i_1, kappa);
        
        c(i*N_x + 1:(i+1)*N_x) = x_i - x_i_1 + dt * (f_i + f_i_1)/2;
    end

% ------------------------------------------------------------------
function J = jacobianstructure(auxdata)  
    [~, ~, ~, ~, ~, N_x, N_u, N_steps, ~] = deal(auxdata{:});

    % Define blocks of full Jacobian
    I = eye(N_x);
    A = [1 1 1 1 1;
         0 1 1 1 1;
         0 1 1 1 1;
         0 0 0 1 0;
         0 0 0 0 1];
    B = [0 0;
         0 0;
         0 0;
         1 0;
         0 1];

    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    J(1:N_x, 1:N_x) = I;
    
    for i = 1:N_steps-1
        J(i*N_x + 1 : (i+1)*N_x, (i-1)*(N_x+N_u)+1 : (i+1)*(N_x+N_u)) = [A, B, A, B];
    end

    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [~, ~, kappa, ~, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    J(1:N_x, 1:N_x) = I;
    
    for i = 1:N_steps-1
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        x_i_1 = x(i*(N_x+N_u) + 1:i*(N_x+N_u) + N_x);
        u_i_1 = x(i*(N_x+N_u) + N_x + 1  :  (i+1)*(N_x+N_u));
        
        A_i = A_curv_kin(x_i, u_i, kappa)/2 * dt + I;
        A_i_1 = A_curv_kin(x_i_1, u_i_1, kappa)/2 * dt - I;
        B_i = B_curv_kin(x_i, u_i, kappa)/2 * dt;
        B_i_1 = B_curv_kin(x_i_1, u_i_1, kappa)/2 * dt;
        
        J(i*N_x + 1 : (i+1)*N_x, (i-1)*(N_x+N_u)+1 : (i+1)*(N_x+N_u)) ...
            = [A_i, B_i, A_i_1, B_i_1];
    end
    
    J = sparse(J);
