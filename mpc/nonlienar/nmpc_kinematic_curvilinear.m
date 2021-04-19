function [x, info] = nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, x_init, info)
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
    Q = [5; 50; 10; 0; 0];
    Q_terminal = Q * 10;
    R = [10, 10];
    
    Q_vec = [repmat([Q(:); R(:)], N_steps - 1, 1); Q_terminal(:); R(:)];
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { x0, x_ref, kappa, Q_bar, N_x, N_u, N_steps, dt };

    % The constraint functions are bounded from below by zero.
    options.lb = repmat([-1e8; -1.0; -1e8; 0; -0.4; -10.0; -0.4], N_steps, 1); % Lower bound on optimization variable
    options.ub = repmat([1e8; 1.0; 1e8; 1e8; 0.4; 10.0; 0.4], N_steps, 1); % Upper bound on optimization variable
    options.cl = zeros(N_x*N_steps, 1); % Lower bound on constraint function
    options.cu = zeros(N_x*N_steps, 1); % Upper bound on constraint function
    
    % Set IPOPT options
    options.ipopt.print_level           = 0;
    options.ipopt.max_iter              = 20;
    options.ipopt.tol                   = 1e-8;
    options.ipopt.hessian_approximation = 'limited-memory';

    % Define callback functions
    funcs.objective         = @objective; % Objective function (Required)
    funcs.gradient          = @gradient; % Gradient of objective (Required)
    
    funcs.constraints       = @constraints; % Constraint function, corresponds to cu/cl (Optional)
    funcs.jacobian          = @jacobian; % Jacobian of the constraints at the current points (Optional)
    funcs.jacobianstructure = @jacobianstructure; %Structure of Jacobian (Optional)

    % Run IPOPT.
    x_init(1:end-(N_x+N_u)*N_steps) = x_init((N_x+N_u)*N_steps+1:end);
    [x, info] = ipopt_auxdata(x_init, funcs, options);    
    
% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [~, x_ref, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x - x_ref(:);
    f = x_error' * Q_bar * x_error;

% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [~, x_ref, ~, Q_bar, ~, ~, ~, ~] = deal(auxdata{:});
    
    x_error = x - x_ref(:);
    g = 2 * Q_bar * x_error;

% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [x0, ~, kappa, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    
    c = zeros(N_x*N_steps, 1);
    
    x_state = x0;

    for i = 1:N_steps
        x_state_next = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        x_control = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        beta = atan(lr/(lr+lf) * tan(x_state(5)));
        k = kappa(x_state(1));
        
        s_d = x_state(4) * cos(x_state(3) + beta) / (1 + x_state(2) * k);
    
        x_state_d = [s_d;
                     x_state(4) * sin(x_state(3) + beta);
                     x_state(4) * sin(beta) / lr - s_d * k;
                     x_control(1);
                     x_control(2)];
        
        c((i-1)*N_x + 1:i*N_x) = x_state + dt * x_state_d - x_state_next;
        x_state = x_state_next;
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
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    J(1:N_x, 1:(N_x+N_u)) = [I, B];    
    
    for i = 2:N_steps
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [I, B];
    end

    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [~, ~, kappa, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = -eye(N_x);

    B = [0 0;
         0 0;
         0 0;
         1 0;
         0 1] * dt;

    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);
    
    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    J(1:N_x, 1:(N_x+N_u)) = [I, B];       
    
    for i = 2:N_steps
        x_state = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        
        % Precompute common expressions
        k = kappa(x_state(1));
        beta = atan(lr_ratio * tan(x_state(5)));
        s_mu_beta = sin(x_state(3) + beta);
        c_mu_beta = cos(x_state(3) + beta);
        beta_d = lr_ratio * sec(x_state(5))^2 / ...
            (1 + (lr_ratio * tan(x_state(5)))^2);
        denom_nk = 1 / (1 - x_state(2) * k); 

        % Compute partial derivatives
        s_n = x_state(4)*c_mu_beta * denom_nk^2 * k;
        s_mu = -x_state(4)*s_mu_beta * denom_nk;
        s_v = c_mu_beta * denom_nk;
        s_delta = -x_state(4)*s_mu_beta * denom_nk * beta_d;

        n_mu = x_state(4)*c_mu_beta;
        n_v = s_mu_beta;
        n_delta = x_state(4)*c_mu_beta * beta_d;

        mu_n = -s_n * k;
        mu_mu = -s_mu * k;
        mu_v = sin(beta)/lr - s_v * k;
        mu_delta = x_state(4)*cos(beta)*beta_d/lr - s_delta * k;

        % Populate matrices
        A = [0   s_n    s_mu    s_v    s_delta;
             0   0      n_mu    n_v    n_delta;
             0   mu_n   mu_mu   mu_v   mu_delta;
             0   0      0       0      0;
             0   0      0       0      0];
        A = A*dt + eye(N_x);
        
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [I, B];
    end
    
    J = sparse(J);
