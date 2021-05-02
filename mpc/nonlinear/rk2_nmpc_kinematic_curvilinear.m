function [x, info] = rk2_nmpc_kinematic_curvilinear(x0, x_ref, kappa, kappa_d, dt, x_init, info)
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
%     options.ipopt.derivative_test_tol   = 1;

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
    
    x_i = x0;

    for i = 1:N_steps
        x_i_1 = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
                
        % Perform RK4 step
        k1 = f_curv_kin(x_i, u_i, kappa);
        k2 = f_curv_kin(x_i + k1*dt / 2, u_i, kappa);
        
        c((i-1)*N_x + 1:i*N_x) = x_i + dt*k2 - x_i_1;
        x_i = x_i_1;
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
    B = [1 1;
         1 1;
         1 1;
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
    [x0, ~, kappa, ~, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    
    u0 = x(N_x + 1 : N_x+N_u);
    B = B_curv_kin(x0, u0, kappa) * dt;
    J(1:N_x, 1:(N_x+N_u)) = [-I, B];       
    
    for i = 2:N_steps
        x_i = x((i-2)*(N_x+N_u) + 1 : (i-2)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : i*(N_x+N_u));
        
        % Calculate RK4 slopes
        k1 = f_curv_kin(x_i, u_i, kappa);
                
        % Calculate RK4 partial derivatives w.r.t. states 
        dfdx1 = A_curv_kin(x_i, u_i, kappa);
        dfdx2 = A_curv_kin(x_i + k1*dt / 2, u_i, kappa);
        
        dkdx1 = dfdx1;
        dkdx2 = dfdx2 * (I + dkdx1*dt / 2);
        
        % Calculate RK4 partial derivatives w.r.t. controls
        dkdu1 = B_curv_kin(x_i, u_i, kappa);
        dkdu2 = B_curv_kin(x_i + k1*dt / 2, u_i, kappa) + dfdx2*dkdu1*dt / 2;
        
        % Calculate final matrices
        A = I + dt * dkdx2;
        B = dt * dkdu2;
        
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [-I, B];
    end
    
    J = sparse(J);
