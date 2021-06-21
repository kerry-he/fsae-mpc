function [x, info, ds, N_steps] = minimum_time_planner(x_spline, y_spline, dl, L)
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

    kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl); 
    
    % Define constants
    N_x = 4;
    N_u = 2;
    N_steps = 1000;
    ds = L / N_steps;

    % Defining cost weights
    Q = [0.01; 0.01; 0; 0];
    Q_terminal = Q * 10;
    R = [1e-6, 0.01];
    
    Q_vec = [repmat([Q(:); R(:)], N_steps - 1, 1); Q_terminal(:); R(:)];
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
        
    % Set up the auxiliary data.
    options.auxdata = { kappa, Q_bar, N_x, N_u, N_steps, ds };

    % The constraint functions are bounded from below by zero.
    options.lb = repmat([-1.0; -inf; 0; -0.4; -10.0; -0.4], N_steps, 1); % Lower bound on optimization variable
    options.ub = repmat([1.0; inf; 50; 0.4; 10.0; 0.4], N_steps, 1); % Upper bound on optimization variable
    options.cl = zeros(N_x*N_steps, 1); % Lower bound on constraint function
    options.cu = zeros(N_x*N_steps, 1); % Upper bound on constraint function
    
    % Set IPOPT options
%     options.ipopt.print_level           = 0;
    options.ipopt.max_iter              = 1000;
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
    x_init = zeros(N_x+N_u, N_steps);
    x_init(3, :) = 20;
    [x, info] = ipopt_auxdata(x_init(:), funcs, options);  
    
% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [kappa, Q_bar, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    f = x' * Q_bar * x;
    
    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    
    for i = 1:N_steps
        s = ds * (i - 1);
        x_state = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        
        beta = atan(lr/(lr+lf) * tan(x_state(4)));
        k = kappa(s);
        
        s_d = x_state(3) * cos(x_state(2) + beta) / (1 - x_state(1) * k);
        
        f = f + ds / s_d;
    end

% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [kappa, Q_bar, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    g = 2 * Q_bar * x;

    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);
    
    for i = 1:N_steps
        x_state = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        
        % Precompute common expressions
        s = ds * (i - 1);
        k = kappa(s);
        beta = atan(lr_ratio * tan(x_state(4)));
        s_mu_beta = sin(x_state(2) + beta);
        c_mu_beta = cos(x_state(2) + beta);
        beta_d = lr_ratio * sec(x_state(4))^2 / ...
            (1 + (lr_ratio * tan(x_state(4)))^2);
        denom_nk = 1 / (1 - x_state(1) * k);         
        s_d = x_state(3) * c_mu_beta * denom_nk;
        
        % Compute partial derivatives
        s_n = x_state(3)*c_mu_beta * denom_nk^2 * k;
        s_mu = -x_state(3)*s_mu_beta * denom_nk;
        s_v = c_mu_beta * denom_nk;
        s_delta = -x_state(3)*s_mu_beta * denom_nk * beta_d;
        
        g((i-1)*(N_x+N_u) + 1 : i*(N_x+N_u)) = g((i-1)*(N_x+N_u) + 1 : i*(N_x+N_u)) - ...
            [s_n;
             s_mu;
             s_v;
             s_delta;
             0;
             0] * ds / s_d^2;
    end

% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [kappa, ~, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    
    c = zeros(N_x*N_steps, 1);

    for i = 1:N_steps
        x_state = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        x_state_next = x(mod(i, N_steps)*(N_x+N_u) + 1:mod(i, N_steps)*(N_x+N_u) + N_x);
        x_control = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        s = ds * (i - 1);
        beta = atan(lr/(lr+lf) * tan(x_state(4)));
        k = kappa(s);
        
        s_d = x_state(3) * cos(x_state(2) + beta) / (1 - x_state(1) * k);
    
        x_state_d = [x_state(3) * sin(x_state(2) + beta);
                     x_state(3) * sin(beta) / lr - s_d * k;
                     x_control(1);
                     x_control(2)] / s_d;
        
        c((i-1)*N_x + 1:i*N_x) = x_state + ds * x_state_d - x_state_next;
    end

% ------------------------------------------------------------------
function J = jacobianstructure(auxdata)  
    [~, ~, N_x, N_u, N_steps, ~] = deal(auxdata{:});

    % Define blocks of full Jacobian
    I = eye(N_x);
    A = ones(4);
    B = [0 0;
         0 0;
         1 0;
         0 1];

    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    
    for i = 1:N_steps
        i_wrap = mod(i, N_steps);
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [A, B];
        J((i-1)*N_x+1 : i*N_x, i_wrap*(N_x+N_u)+1 : i_wrap*(N_x+N_u)+N_x) = I;
    end

    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [kappa, ~, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    % Define vehicle constants
    lr = 0.6183;
    lf = 0.8672;
    lr_ratio = lr / (lr + lf);

    I = -eye(N_x);
    
    % Fill out Jacobian
    J = zeros(N_x*N_steps, (N_x+N_u)*N_steps);
    
    for i = 1:N_steps
        s = ds * i;
        x_state = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        x_control = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        % Precompute common expressions
        k = kappa(s);
        beta = atan(lr_ratio * tan(x_state(4)));
        s_mu_beta = sin(x_state(2) + beta);
        c_mu_beta = cos(x_state(2) + beta);
        beta_d = lr_ratio * sec(x_state(4))^2 / ...
            (1 + (lr_ratio * tan(x_state(4)))^2);
        denom_nk = 1 / (1 - x_state(1) * k); 

        % Compute partial derivatives
        s_d = x_state(3) * c_mu_beta * denom_nk;
        f = [x_state(3) * sin(x_state(2) + beta);
             x_state(3) * sin(beta) / lr - s_d * k;
             x_control(1);
             x_control(2)];        
        
        sd_x = [x_state(3)*c_mu_beta * denom_nk^2 * k;
                -x_state(3)*s_mu_beta * denom_nk;
                c_mu_beta * denom_nk;
                -x_state(3)*s_mu_beta * denom_nk * beta_d];        
        
        n_mu = x_state(3)*c_mu_beta;
        n_v = s_mu_beta;
        n_delta = x_state(3)*c_mu_beta * beta_d;

        mu_n = -sd_x(1) * k;
        mu_mu = -sd_x(2) * k;
        mu_v = sin(beta)/lr - sd_x(3) * k;
        mu_delta = x_state(3)*cos(beta)*beta_d/lr - sd_x(4) * k;

        % Populate matrices
        A = [0      n_mu    n_v    n_delta;
             mu_n   mu_mu   mu_v   mu_delta;
             0      0       0      0;
             0      0       0      0] / s_d ...
            - f * sd_x' / s_d^2 ;
        A = A*ds + eye(N_x);

        B = [0 0;
             0 0;
             1 0;
             0 1] * ds / s_d;        
        
        i_wrap = mod(i, N_steps);
        J((i-1)*N_x+1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [A, B];
        J((i-1)*N_x+1 : i*N_x, i_wrap*(N_x+N_u)+1 : i_wrap*(N_x+N_u)+N_x) = I;
    end
    
    J = sparse(J);
