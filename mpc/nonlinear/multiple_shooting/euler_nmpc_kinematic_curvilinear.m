function [x, info] = euler_nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, x_init, info)
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
    N_u = 2;
    
    x_ref = [x_ref; zeros(N_u, N_steps)];

    % Defining cost weights
    Q = [5; 250; 2000; 0; 0];
    Q_terminal = Q * 10;
    R = [10, 10];
    
    Q_vec = [repmat([Q(:); R(:)], N_steps - 1, 1); Q_terminal(:); R(:)];
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { x0, x_ref, kappa, Q_bar, N_x, N_u, N_steps, dt };

    % The constraint functions are bounded from below by zero.
    options.lb = [repmat([-inf; -inf; -inf; 0; -0.4; -10.0; -0.4], N_steps, 1); 0]; % Lower bound on optimization variable
    options.ub = [repmat([inf; inf; inf; inf; 0.4; 10.0; 0.4], N_steps, 1); inf]; % Upper bound on optimization variable
    options.cl = [zeros(N_x*N_steps, 1); repmat([-inf; -0.75], N_steps, 1); repmat([-inf; -5.0], N_steps, 1)]; % Lower bound on constraint function
    options.cu = [zeros(N_x*N_steps, 1); repmat([0.75; inf], N_steps, 1); repmat([5.0; inf], N_steps, 1)]; % Upper bound on constraint function
    
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
    x_init(1:(N_x+N_u)*(N_steps-1)) = x_init(N_x+N_u+1:(N_x+N_u)*N_steps);
    x_init(end-(N_x+N_u) : end-N_u-1) = x_init(end-(N_x+N_u) : end-N_u-1)...
        + dt*f_curv_kin(x_init(end-(N_x+N_u) : end-N_u-1), x_init(end-N_u : end-1), kappa);
    x_init(end) = 0;
    [x, info] = ipopt_auxdata(x_init(:), funcs, options);  
    
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
    c = zeros((N_x + 2 + 2)*N_steps, 1);
    
    x_i = x0;

    for i = 1:N_steps
        x_i_1 = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        % Calculate dynamic model
        f = f_curv_kin(x_i, u_i, kappa);
        
        c((i-1)*N_x + 1:i*N_x) = x_i + dt * f - x_i_1;
        x_i = x_i_1;
        
        
        % Soft constraints
        c(N_x*N_steps + 1 + 2*(i-1)) = x_i(2) - x(end);
        c(N_x*N_steps + 2 + 2*(i-1)) = x_i(2) + x(end);
        
        
        % Friction constraints
        lr = 0.6183;
        lf = 0.8672;        
        c(N_x*N_steps + 2*N_steps + 1 + 2*(i-1)) = x_i(4)^2 * x_i(5) / (lr + lf) - x(end);
        c(N_x*N_steps + 2*N_steps + 2 + 2*(i-1)) = x_i(4)^2 * x_i(5) / (lr + lf) + x(end);    
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
    J = zeros((N_x + 2 + 2)*N_steps, (N_x+N_u)*N_steps + 1);
    J(1:N_x, 1:(N_x+N_u)) = [I, B];    
    
    for i = 2:N_steps
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [I, B];
    end
    
    
    % Slack constraints
    for i = 1:N_steps
        J(N_x*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 2) = 1;
        J(N_x*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 2) = 1;        
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*N_steps, end) = ones(N_steps*2, 1);
    
    
    %  Friction constraints
    for i = 1:N_steps    
        J(N_x*N_steps + 2*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 4) = 1;  
        J(N_x*N_steps + 2*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 5) = 1;  
        
        J(N_x*N_steps + 2*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 4) = 1;  
        J(N_x*N_steps + 2*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 5) = 1;          
    end        
    
    J(N_x*N_steps + 2*N_steps + 1 : end, end) = ones(N_steps*2, 1);        

    
    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [x0, ~, kappa, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros((N_x + 2 + 2)*N_steps, (N_x+N_u)*N_steps + 1);
    
    u0 = x(N_x + 1 : N_x+N_u);
    B = B_curv_kin(x0, u0, kappa) * dt;
    J(1:N_x, 1:(N_x+N_u)) = [-I, B];       
    
    for i = 2:N_steps
        x_i = x((i-2)*(N_x+N_u) + 1 : (i-2)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : i*(N_x+N_u));
        
        A = A_curv_kin(x_i, u_i, kappa) * dt + I;
        B = B_curv_kin(x_i, u_i, kappa) * dt;
        
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [-I, B];
    end
    
    
    % Slack constraints
    for i = 1:N_steps
        J(N_x*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 2) = 1;
        J(N_x*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 2) = 1;        
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*N_steps, end) = repmat([-1; 1], N_steps, 1);    
    
    
    %  Friction constraints
    for i = 1:N_steps
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        
        lr = 0.6183;
        lf = 0.8672;        
        J(N_x*N_steps + 2*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 4) = 2*x_i(4) * x_i(5) / (lr + lf);  
        J(N_x*N_steps + 2*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 5) = x_i(4)^2 / (lr + lf);  
        
        J(N_x*N_steps + 2*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 4) = 2*x_i(4) * x_i(5) / (lr + lf);  
        J(N_x*N_steps + 2*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 5) = x_i(4)^2 / (lr + lf);  
    end      
    
    J(N_x*N_steps + 2*N_steps + 1 : end, end) = repmat([-1; 1], N_steps, 1);  
       
    
    J = sparse(J);
