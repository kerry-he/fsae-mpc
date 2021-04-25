function [x, info] = nmpc_kinematic_curvilinear(x0, x_ref, kappa, kappa_d, dt, x_init, info)
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
    options.ipopt.tol                   = 1e-8; % OR 1e-5
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
    x_init(1:end-(N_x+N_u)*N_steps) = x_init((N_x+N_u)*N_steps+1:end);
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
    
    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    
    c = zeros(N_x*N_steps, 1);
    
    x_state = x0;

    rk_factor = [1 2 2 1];
    for i = 1:N_steps
        x_state_next = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        x_control = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        x_rk = x_state;
        f = zeros(N_x, 1);
        for j = 1:4
            beta = atan(lr/(lr+lf) * tan(x_rk(5)));
            k = kappa(x_rk(1));

            s_d = x_rk(4) * cos(x_rk(3) + beta) / (1 - x_rk(2) * k);

            x_state_d = [s_d;
                         x_rk(4) * sin(x_rk(3) + beta);
                         x_rk(4) * sin(beta) / lr - s_d * k;
                         x_control(1);
                         x_control(2)];
                     
            x_rk = x_state + x_state_d*dt / rk_factor(mod(j, 4)+1);
            f = f + x_state_d*rk_factor(j);
        end
        
        c((i-1)*N_x + 1:i*N_x) = x_state + dt * f / 6 - x_state_next;
        x_state = x_state_next;
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
    J(1:N_x, 1:(N_x+N_u)) = [I, B];    
    
    for i = 2:N_steps
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [I, B];
    end

    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [~, ~, kappa, kappa_d, ~, N_x, N_u, N_steps, dt] = deal(auxdata{:});
    
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
    
    rk_factor = [1 2 2 1];
    for i = 2:N_steps
        x_state = x((i-2)*(N_x+N_u) + 1:(i-2)*(N_x+N_u) + N_x);
        x_control = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        x_k = x_state;
        A = zeros(N_x, N_x);
        dk_dx = zeros(N_x, N_x, 4);
        for j = 1:4
            % Precompute common expressions
            k = kappa(x_k(1));
            k_d = kappa_d(x_k(1));
            beta = atan(lr_ratio * tan(x_k(5)));
            s_mu_beta = sin(x_k(3) + beta);
            c_mu_beta = cos(x_k(3) + beta);
            beta_d = lr_ratio * sec(x_k(5))^2 / ...
                (1 + (lr_ratio * tan(x_k(5)))^2);
            denom_nk = 1 / (1 - x_k(2) * k); 

            % Compute partial derivatives
            s_s = x_k(4)*c_mu_beta * denom_nk^2 * k_d * x_k(2);
            s_n = x_k(4)*c_mu_beta * denom_nk^2 * k;
            s_mu = -x_k(4)*s_mu_beta * denom_nk;
            s_v = c_mu_beta * denom_nk;
            s_delta = -x_k(4)*s_mu_beta * denom_nk * beta_d;

            n_mu = x_k(4)*c_mu_beta;
            n_v = s_mu_beta;
            n_delta = x_k(4)*c_mu_beta * beta_d;

            mu_s = -x_k(4)*c_mu_beta * denom_nk * k_d - s_s * k;
            mu_n = -s_n * k;
            mu_mu = -s_mu * k;
            mu_v = sin(beta)/lr - s_v * k;
            mu_delta = x_k(4)*cos(beta)*beta_d/lr - s_delta * k;

            % Populate matrices
            A_k = [0      s_n    s_mu    s_v    s_delta;
                   0      0      n_mu    n_v    n_delta;
                   0      mu_n   mu_mu   mu_v   mu_delta;
                   0      0      0       0      0;
                   0      0      0       0      0];
               
            switch j
                case 1
                    dk_dx(:, :, j) = A_k;
                case 2
                    dk_dx(:, :, j) = A_k * (eye(N_x) + dk_dx(:, :, j - 1)*dt / 2);
                case 3
                    dk_dx(:, :, j) = A_k * (eye(N_x) + dk_dx(:, :, j - 1)*dt / 2);
                case 4
                    dk_dx(:, :, j) = A_k * (eye(N_x) + dk_dx(:, :, j - 1)*dt);
            end         
               
            A = A + A_k*rk_factor(j);

            f = [x_k(4)*c_mu_beta*denom_nk;
                 x_k(4)*s_mu_beta;
                 x_k(4)*sin(beta)/lr - x_k(4)*c_mu_beta*denom_nk*k;
                 x_control(1);
                 x_control(2)];

            x_k = x_state + f*dt / rk_factor(mod(j, 4)+1);
        end

        % Populate matrices
        A = A*dt / 6 + eye(N_x);
        
        J((i-1)*N_x+1 : i*N_x, (i-2)*(N_x+N_u)+1 : (i-2)*(N_x+N_u)+N_x) = A;
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [I, B];
    end
    
    J = sparse(J);
