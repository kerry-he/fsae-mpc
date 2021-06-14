function [x, slack, info] = trapezoidal_nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, x_init, info)
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
    Q = [5; 250; 2000; 0; 0];
    Q_terminal = 10;
    R = [10; 10];
    
    Q_factor = [1/2, ones(1, N_steps-3), (1+Q_terminal)/2, Q_terminal/2];
    R_factor = [1/2, ones(1, N_steps-2), 1/2];
    Q_factor = repmat(Q_factor, N_x, 1);
    R_factor = repmat(R_factor, N_u, 1);
    
    Q = Q .* Q_factor;
    R = R .* R_factor;
    
    Q_vec = [Q; R]; Q_vec = Q_vec(:);
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { x0, x_ref, kappa, Q_bar, N_x, N_u, N_steps, dt };

    % The constraint functions are bounded from below by zero.
    options.lb = [repmat([-inf; -inf; -inf; 0; -0.4; -10.0; -0.4], N_steps, 1); 0]; % Lower bound on optimization variable
    options.ub = [repmat([inf; inf; inf; inf; 0.4; 10.0; 0.4], N_steps, 1); inf]; % Upper bound on optimization variable
    options.cl = [zeros(N_x*N_steps, 1); repmat([-inf; -0.75], N_steps-1, 1); repmat([-inf; -5.0], N_steps-1, 1)]; % Lower bound on constraint function
    options.cu = [zeros(N_x*N_steps, 1); repmat([0.75; inf], N_steps-1, 1); repmat([5.0; inf], N_steps-1, 1)]; % Upper bound on constraint function
    
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
    x_init(end-(N_x+N_u)+1 : end-N_u) = x_init(end-(N_x+N_u)+1 : end-N_u)...
        + dt*f_curv_kin(x_init(end-(N_x+N_u)+1 : end-N_u), x_init(end-N_u+1 : end), kappa);
    x_init = [x_init; 0];
    [x, info] = ipopt_auxdata(x_init(:), funcs, options);  
    
    slack = x(end);
    x = x(1:end-1);
    
% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [~, x_ref, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x(1:end-1) - x_ref(:);
    f = x_error' * Q_bar * x_error + x(end)*1e8;

% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [~, x_ref, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x(1:end-1) - x_ref(:);
    g = [2 * Q_bar * x_error; 1e8];

% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [x0, ~, kappa, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    % Preallocate
     c = zeros((N_x + 2 + 2)*N_steps - 4, 1);
    
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
    
    for i = 2:N_steps
        % Soft constraints
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        c(N_x*N_steps + 1 + 2*(i-2)) = x_i(2) - x(end);
        c(N_x*N_steps + 2 + 2*(i-2)) = x_i(2) + x(end);                   
        
        % Friction constraints
        lr = 0.6183;
        lf = 0.8672;        
        c(N_x*N_steps + 2*N_steps - 1 + 2*(i-2)) = x_i(4)^2 * x_i(5) / (lr + lf) - x(end);
        c(N_x*N_steps + 2*N_steps + 2*(i-2)) = x_i(4)^2 * x_i(5) / (lr + lf) + x(end);
    end

% ------------------------------------------------------------------
function J = jacobianstructure(auxdata)  
    [~, ~, ~, ~, N_x, N_u, N_steps, ~] = deal(auxdata{:});

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
    J = zeros((N_x + 2 + 2)*N_steps - 4, (N_x+N_u)*N_steps + 1);
    J(1:N_x, 1:N_x) = I;
    
    for i = 1:N_steps-1
        J(i*N_x + 1 : (i+1)*N_x, (i-1)*(N_x+N_u)+1 : (i+1)*(N_x+N_u)) = [A, B, A, B];
    end
    
    % Slack constraints
    for i = 2:N_steps
        J(N_x*N_steps + 1 + 2*(i-2), (i-1)*(N_x+N_u) + 2) = 1;
        J(N_x*N_steps + 2 + 2*(i-2), (i-1)*(N_x+N_u) + 2) = 1;  
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*(N_steps-1), end) = ones((N_steps-1)*2, 1);
    
    
    %  Friction constraints
    for i = 2:N_steps     
        J(N_x*N_steps + 2*N_steps - 1 + 2*(i-2), (i-1)*(N_x+N_u) + 4) = 1;  
        J(N_x*N_steps + 2*N_steps - 1 + 2*(i-2), (i-1)*(N_x+N_u) + 5) = 1;  
        
        J(N_x*N_steps + 2*N_steps + 2*(i-2), (i-1)*(N_x+N_u) + 4) = 1;  
        J(N_x*N_steps + 2*N_steps + 2*(i-2), (i-1)*(N_x+N_u) + 5) = 1;          
    end        
    
    J(N_x*N_steps + 2*(N_steps-1) + 1 : end, end) = ones((N_steps-1)*2, 1);  
    
    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [~, ~, kappa, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros((N_x + 2 + 2)*N_steps - 4, (N_x+N_u)*N_steps + 1);
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
    
    % Slack constraints
    for i = 2:N_steps
        J(N_x*N_steps + 1 + 2*(i-2), (i-1)*(N_x+N_u) + 2) = 1;
        J(N_x*N_steps + 2 + 2*(i-2), (i-1)*(N_x+N_u) + 2) = 1;        
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*(N_steps-1), end) = repmat([-1; 1], N_steps-1, 1);      

    
    %  Friction constraints
    for i = 2:N_steps
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        
        lr = 0.6183;
        lf = 0.8672;        
        J(N_x*N_steps + 2*N_steps - 1 + 2*(i-2), (i-1)*(N_x+N_u) + 4) = 2*x_i(4) * x_i(5) / (lr + lf);  
        J(N_x*N_steps + 2*N_steps - 1 + 2*(i-2), (i-1)*(N_x+N_u) + 5) = x_i(4)^2 / (lr + lf);  
        
        J(N_x*N_steps + 2*N_steps + 2*(i-2), (i-1)*(N_x+N_u) + 4) = 2*x_i(4) * x_i(5) / (lr + lf);  
        J(N_x*N_steps + 2*N_steps + 2*(i-2), (i-1)*(N_x+N_u) + 5) = x_i(4)^2 / (lr + lf);  
    end      
    
    J(N_x*N_steps + 2*(N_steps-1) + 1 : end, end) = repmat([-1; 1], N_steps-1, 1);  
    
    J = sparse(J);
