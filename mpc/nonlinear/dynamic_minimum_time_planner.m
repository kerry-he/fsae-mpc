function [x, t, info, ds, N_steps, slack] = dynamic_minimum_time_planner(x_spline, y_spline, dl, L)
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
    N_x = 6;
    N_u = 2;
    N_steps = 500;
    ds = L / N_steps;

    % Defining cost weights
    Q = [0.01; 0.01; 0; 0; 0; 0];
    Q_terminal = 1;
    R = [1e-6; 0.01];
    
    Q_factor = [1/2, ones(1, N_steps-3), (1+Q_terminal)/2, Q_terminal/2];
    R_factor = [1/2, ones(1, N_steps-2), 1/2];
    Q_factor = repmat(Q_factor, N_x, 1);
    R_factor = repmat(R_factor, N_u, 1);
    
    Q = Q .* Q_factor;
    R = R .* R_factor;
    
    Q_vec = [Q; R]; Q_vec = Q_vec(:);
    Q_bar = spdiags(Q_vec, 0, length(Q_vec), length(Q_vec));
    
    % Set up the auxiliary data.
    options.auxdata = { kappa, Q_bar, N_x, N_u, N_steps, ds };

    % The constraint functions are bounded from below by zero.
    options.lb = [repmat([-inf; -inf; 0; -inf; -inf; -0.4; -10.0; -0.4], N_steps, 1); 0; 0]; % Lower bound on optimization variable
    options.ub = [repmat([inf; inf; inf; inf; inf; 0.4; 10.0; 0.4], N_steps, 1); inf; inf]; % Upper bound on optimization variable
    options.cl = [zeros(N_x*N_steps, 1); repmat([-inf; -0.5], N_steps, 1); -inf*ones(N_steps, 1)]; % Lower bound on constraint function
    options.cu = [zeros(N_x*N_steps, 1); repmat([0.5; inf], N_steps, 1); 0.75*ones(N_steps, 1)]; % Upper bound on constraint function
    
    % Set IPOPT options
%     options.ipopt.print_level           = 0;
    options.ipopt.max_iter              = 1000;
    options.ipopt.tol                   = 1e-5; % OR 1e-5
    options.ipopt.hessian_approximation = 'limited-memory';
%     options.ipopt.derivative_test       = 'first-order';
%     options.ipopt.derivative_test_tol   = 1e-0;
    
    % Define callback functions
    funcs.objective         = @objective; % Objective function (Required)
    funcs.gradient          = @gradient; % Gradient of objective (Required)
    
    funcs.constraints       = @constraints; % Constraint function, corresponds to cu/cl (Optional)
    funcs.jacobian          = @jacobian; % Jacobian of the constraints at the current points (Optional)
    funcs.jacobianstructure = @jacobianstructure; %Structure of Jacobian (Optional)

    % Run IPOPT.
    x_init = zeros(N_x+N_u, N_steps);
    x_init(3, :) = 10;
    x_init = [x_init(:); 0; 0];
    [x, info] = ipopt_auxdata(x_init, funcs, options);  
    
    slack = x(end-1:end);
    x = x(1:end-2);
    
    % Calculate time vector
    t = zeros(N_steps, 1);
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        state_d = f_curv_dyn([s; x_i], u_i, kappa);
        s_d = state_d(1);
        
        t(i) = ds / s_d;
    end    
    
% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [kappa, Q_bar, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    % Reference cost
    x_error = x(1:end-2);
    f = x_error' * Q_bar * x_error + x(end-1)*1e5 + x(end)*1e8;
    
    % Lap time cost    
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        state_d = f_curv_dyn([s; x_i], u_i, kappa);
        s_d = state_d(1);
        
        f = f + ds / s_d;
    end    

% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [kappa, Q_bar, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    % Reference cost
    x_error = x(1:end-2);
    g = [2 * Q_bar * x_error; 1e5; 1e8];
    
    % Lap time cost  
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        
        % Precompute common expressions
        f = f_curv_dyn([s; x_i], u_i, kappa);
        s_d = f(1);
        
        A = A_curv_dyn([s; x_i], u_i, kappa);
        s_partial = A(1, 2:end);
        
        g((i-1)*(N_x+N_u) + 1 : i*(N_x+N_u)) = g((i-1)*(N_x+N_u) + 1 : i*(N_x+N_u)) - ...
            [s_partial'; 0; 0] * ds / s_d^2;
    end    

% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [kappa, ~, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    % Preallocate
    c = zeros((N_x + 2 + 1)*N_steps, 1);
        
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : i*(N_x+N_u));
        x_i_1 = x(mod(i, N_steps)*(N_x+N_u) + 1 : mod(i, N_steps)*(N_x+N_u) + N_x);
        u_i_1 = x(mod(i, N_steps)*(N_x+N_u) + N_x + 1 : mod(i, N_steps)*(N_x+N_u) + N_x + N_u);
        
        % Calculate dynamic model
        f_i = f_curv_dyn([s; x_i], u_i, kappa);
        f_i_1 = f_curv_dyn([s+ds; x_i_1], u_i_1, kappa);
        
        c((i-1)*N_x + 1:i*N_x) = x_i - x_i_1 + ds * (f_i(2:end)/f_i(1) + f_i_1(2:end)/f_i_1(1))/2;          
    end
    
    for i = 1:N_steps
        % Soft constraints
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        c(N_x*N_steps + 1 + 2*(i-1)) = x_i(1) - x(end);
        c(N_x*N_steps + 2 + 2*(i-1)) = x_i(1) + x(end);                   
    end
    
    for i = 1:N_steps
        % Friction constraints
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        [~, Fcr] = f_curv_dyn([s; x_i], u_i, kappa);
        
        ac_max = 9.163;
        al_max = 10.0;        
        c(N_x*N_steps + 2*N_steps + i) = (Fcr / (280*ac_max))^2 + (u_i(1) / al_max)^2 - x(end-1);                
    end    
    
% ------------------------------------------------------------------
function J = jacobianstructure(auxdata)  
    [~, ~, N_x, N_u, N_steps, ~] = deal(auxdata{:});

    % Define blocks of full Jacobian
    A = ones(6);
     
    B = [0 0;
         0 0;
         1 0;
         0 0;
         0 0;
         0 1];

    % Fill out Jacobian
    J = zeros((N_x + 2 + 1)*N_steps, (N_x+N_u)*N_steps + 2);
    
    for i = 1:N_steps
        J((i-1)*N_x + 1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [A, B];
        J((i-1)*N_x + 1 : i*N_x, mod(i, N_steps)*(N_x+N_u)+1 : mod(i, N_steps)*(N_x+N_u) + N_x + N_u) = [A, B];
    end
    
    % Slack constraints
    for i = 1:N_steps
        J(N_x*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 1) = 1;
        J(N_x*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 1) = 1;  
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*N_steps, end) = ones(N_steps*2, 1);
    
    
    % Friction constraints
    for i = 1:N_steps    
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 3) = 1;  
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 4) = 1; 
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 5) = 1;  
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 7) = 1;   
    end        
    
    J(N_x*N_steps + 2*N_steps + 1 : end, end-1) = ones(N_steps, 1);   
    
    
    J = sparse(J);

% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    [kappa, ~, N_x, N_u, N_steps, ds] = deal(auxdata{:});
    
    I = eye(N_x);

    % Fill out Jacobian
    J = zeros((N_x + 2 + 1)*N_steps, (N_x+N_u)*N_steps + 2);
    
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1 : (i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1 : i*(N_x+N_u));
        x_i_1 = x(mod(i, N_steps)*(N_x+N_u) + 1 : mod(i, N_steps)*(N_x+N_u) + N_x);
        u_i_1 = x(mod(i, N_steps)*(N_x+N_u) + N_x + 1 : mod(i, N_steps)*(N_x+N_u) + N_x + N_u);
        
        % Calculate dynamic model
        f_i = f_curv_dyn([s; x_i], u_i, kappa);
        f_i_1 = f_curv_dyn([s+ds; x_i_1], u_i_1, kappa);
        
        A_i = A_curv_dyn([s; x_i], u_i, kappa);
        A_i_1 = A_curv_dyn([s+ds; x_i_1], u_i_1, kappa);
        B_i = B_curv_dyn([s; x_i], u_i, kappa);
        B_i_1 = B_curv_dyn([s+ds; x_i_1], u_i_1, kappa);
        
        A_i = (A_i(2:end, 2:end)/f_i(1) - f_i(2:end)*A_i(1, 2:end)/f_i(1)^2)/2 * ds + I;
        A_i_1 = (A_i_1(2:end, 2:end)/f_i_1(1) - f_i_1(2:end)*A_i_1(1, 2:end)/f_i_1(1)^2)/2 * ds - I;
        B_i = B_i(2:end, :)/f_i(1) /2 * ds;
        B_i_1 = B_i_1(2:end, :)/f_i_1(1) /2 * ds;
        
        J((i-1)*N_x + 1 : i*N_x, (i-1)*(N_x+N_u)+1 : i*(N_x+N_u)) = [A_i, B_i];
        J((i-1)*N_x + 1 : i*N_x, mod(i, N_steps)*(N_x+N_u)+1 : mod(i, N_steps)*(N_x+N_u) + N_x + N_u) = [A_i_1, B_i_1];
    end
    
    
    % Slack constraints
    for i = 1:N_steps
        J(N_x*N_steps + 1 + 2*(i-1), (i-1)*(N_x+N_u) + 1) = 1;
        J(N_x*N_steps + 2 + 2*(i-1), (i-1)*(N_x+N_u) + 1) = 1;       
    end

    J(N_x*N_steps + 1 : N_x*N_steps + 2*N_steps, end) = repmat([-1; 1], N_steps, 1);      
    
    
    % Friction constraints
    for i = 1:N_steps
        s = ds * (i - 1);
        x_i = x((i-1)*(N_x+N_u) + 1:(i-1)*(N_x+N_u) + N_x);
        u_i = x((i-1)*(N_x+N_u) + N_x + 1  :  i*(N_x+N_u));
        [~, Fcr, Fcr_d, vr, denom_vr2, x_d_hat, x_d_hat_d] = A_curv_dyn([s; x_i], u_i, kappa);
        
        % Friction constraints
        ac_max = 9.163;
        al_max = 10.0; 
        lr = 0.6183;    
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 3) = 2*Fcr*Fcr_d*denom_vr2*vr*x_d_hat_d/x_d_hat / (280*ac_max)^2;  
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 4) = -2*Fcr*Fcr_d*denom_vr2/x_d_hat / (280*ac_max)^2; 
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 5) = 2*Fcr*Fcr_d*denom_vr2*lr/x_d_hat / (280*ac_max)^2;  
        J(N_x*N_steps + 2*N_steps + i, (i-1)*(N_x+N_u) + 7) = 2*u_i(1) / al_max^2; 
    end      
    
    J(N_x*N_steps + 2*N_steps + 1 : end, end-1) = -ones(N_steps, 1);    

    
    J = sparse(J);
