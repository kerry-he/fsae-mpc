function [x, info] = hs_nmpc_kinematic_curvilinear(x0, x_ref, kappa, kappa_d, dt, x_init, info)
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
    x_ref_interp = zeros(N_x+N_u, N_steps*2-1);
    for i = 1:length(x_ref)-1
        x_ref_interp(:, (i-1)*2 + 1) = x_ref(:, i);
        x_ref_interp(:, (i-1)*2 + 2) = (x_ref(:, i) + x_ref(:, i+1)) / 2;
        x_ref_interp(:, (i-1)*2 + 3) = x_ref(:, i+1);
    end
    x_ref = x_ref_interp;

    % Defining cost weights
    Q = [5; 500; 2000; 0; 0];
    Q_terminal = Q * 10;
    R = [10, 10];
    
    Q_vec = [repmat([Q(:); R(:); Q(:)*4; R(:)*4], N_steps - 1, 1); Q_terminal(:); R(:)];
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { x0, x_ref, kappa, kappa_d, Q_bar, N_x, N_u, N_steps, dt };

    % The constraint functions are bounded from below by zero.
    options.lb = repmat([-inf; -1.0; -inf; 0; -0.4; -10.0; -0.4], N_steps*2-1, 1); % Lower bound on optimization variable
    options.ub = repmat([inf; 1.0; inf; inf; 0.4; 10.0; 0.4], N_steps*2-1, 1); % Upper bound on optimization variable
    options.cl = zeros(N_x*(2*N_steps-1), 1); % Lower bound on constraint function
    options.cu = zeros(N_x*(2*N_steps-1), 1); % Upper bound on constraint function
    
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
    
    for i = 1:2:(N_steps-1)*2
        % Extract required states and controls from state vector
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : (i-1)*(N_x+N_u) + N_x + N_u);
        
        x_i_c = x(i*(N_x+N_u) + 1 : i*(N_x+N_u) + N_x);
        u_i_c = x(i*(N_x+N_u) + N_x + 1 : i*(N_x+N_u) + N_x + N_u);
        
        x_i_1 = x((i+1)*(N_x+N_u) + 1 : (i+1)*(N_x+N_u) + N_x);
        u_i_1 = x((i+1)*(N_x+N_u) + N_x + 1 : (i+2)*(N_x+N_u));
        
        % Calculate dynamic model
        f_i = f_curv_kin(x_i, u_i, kappa);
        f_i_c = f_curv_kin(x_i_c, u_i_c, kappa);
        f_i_1 = f_curv_kin(x_i_1, u_i_1, kappa);
        
        c(i*N_x + 1:(i+1)*N_x) = (x_i + x_i_1) / 2 + (f_i - f_i_1)*dt / 8 - x_i_c;
        c((i+1)*N_x + 1:(i+2)*N_x) = x_i - x_i_1 + (f_i + 4*f_i_c + f_i_1)*dt / 6;
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
    J = zeros(N_x*(2*(N_steps-1) + 1), (N_x+N_u)*(2*N_steps-1));
    J(1:N_x, 1:N_x) = I;
    
    for i = 1:2:(N_steps-1)*2
        J(i*N_x + 1 : (i+1)*N_x, (i-1)*(N_x+N_u)+1 : (i+2)*(N_x+N_u)) ...
            = [A, B, I, zeros(N_x, N_u), A, B];
        
        J((i+1)*N_x + 1 : (i+2)*N_x, (i-1)*(N_x+N_u)+1 : (i+2)*(N_x+N_u))...
            = [A, B, A, B, A, B];
    end

    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [~, ~, kappa, ~, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros(N_x*(2*(N_steps-1) + 1), (N_x+N_u)*(2*N_steps-1));
    J(1:N_x, 1:N_x) = I;
    
    for i = 1:2:(N_steps-1)*2
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : (i-1)*(N_x+N_u) + N_x + N_u);
        
        x_i_c = x(i*(N_x+N_u) + 1 : i*(N_x+N_u) + N_x);
        u_i_c = x(i*(N_x+N_u) + N_x + 1 : i*(N_x+N_u) + N_x + N_u);
        
        x_i_1 = x((i+1)*(N_x+N_u) + 1 : (i+1)*(N_x+N_u) + N_x);
        u_i_1 = x((i+1)*(N_x+N_u) + N_x + 1 : (i+2)*(N_x+N_u));
        
        A_i = A_curv_kin(x_i, u_i, kappa);
        A_i_c = A_curv_kin(x_i_c, u_i_c, kappa);
        A_i_1 = A_curv_kin(x_i_1, u_i_1, kappa);
        
        B_i = B_curv_kin(x_i, u_i, kappa);
        B_i_c = B_curv_kin(x_i_c, u_i_c, kappa);
        B_i_1 = B_curv_kin(x_i_1, u_i_1, kappa);
        
        
        J(i*N_x + 1 : (i+1)*N_x, (i-1)*(N_x+N_u)+1 : (i+2)*(N_x+N_u)) ...
            = [I/2 + A_i*dt/8, B_i*dt/8, -I, zeros(N_x, N_u), I/2 - A_i_1*dt/8, -B_i_1*dt/8];
        
        J((i+1)*N_x + 1 : (i+2)*N_x, (i-1)*(N_x+N_u)+1 : (i+2)*(N_x+N_u))...
            = [I + A_i*dt/6, B_i*dt/6, A_i_c*dt*2/3, B_i_c*dt*2/3, -I + A_i_1*dt/6, B_i_1*dt/6];
    end
    
    J = sparse(J);
