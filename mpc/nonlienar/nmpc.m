function w = nmpc(A, y, lambda)

    % Get the number of examples (n) and the number of regression
    % coefficients (m).
    [n m] = size(A);

    % The starting point.
    x0 = { zeros(m,1) 
         ones(m,1) };         

    % Set up the auxiliary data.
    options.auxdata = { n m A y lambda };

    % The constraint functions are bounded from below by zero.
    options.lb = []; % Lower bound on optimization variable
    options.ub = []; % Upper bound on optimization variable
    options.cl = zeros(2*m,1); % Lower bound on constraint function
    options.cu = repmat(Inf,2*m,1); % Upper bound on constraint function
    
    % Set IPOPT options
    options.ipopt.print_level      = 0;
    options.ipopt.max_iter         = 100;
    options.ipopt.tol              = 1e-8;
    options.ipopt.derivative_test  = 'first-order'; % Debugging

    % Define callback functions
    funcs.objective         = @objective; % Objective function (Required)
    funcs.gradient          = @gradient; % Gradient of objective (Required)
    funcs.hessian           = @hessian; % Hessian of the Lagrangian at the current point (Required unless using quasi-Newton approximation)
    funcs.hessianstructure  = @hessianstructure; % Structure of Hessian (Required unless using quasi-Newton approximation)
    
    funcs.constraints       = @constraints; % Constraint function, corresponds to cu/cl (Optional)
    funcs.jacobian          = @jacobian; % Jacobian of the constraints at the current points (Optional)
    funcs.jacobianstructure = @jacobianstructure; %Structure of Jacobian (Optional)


    % Run IPOPT.
    [x info] = ipopt_auxdata(x0, funcs, options);
    w        = x{1};

% ------------------------------------------------------------------
function f = objective(x, auxdata)
    [n m A y lambda] = deal(auxdata{:});
    [w u] = deal(x{:});
    f     = norm(y - A*w)^2/2 + lambda*sum(u);
  
% ------------------------------------------------------------------
function c = constraints(x, auxdata)
    [w u] = deal(x{:});
    c     = [ w + u; u - w ];
  
% ------------------------------------------------------------------
function g = gradient(x, auxdata)
    [n m A y lambda] = deal(auxdata{:});
    w = x{1};
    g = { -A'*(y - A*w) 
        repmat(lambda,m,1) };
  
% ------------------------------------------------------------------
function J = jacobianstructure(auxdata)  
    m = auxdata{2};
    I = speye(m);
    J = [ I I
        I I ];
  
% ------------------------------------------------------------------
function J = jacobian(x, auxdata)  
    m = auxdata{2};
    I = speye(m);
    J = [  I  I
        -I  I ];
  
% ------------------------------------------------------------------
function H = hessianstructure(auxdata)
    m = auxdata{2};
    H = [ tril(ones(m))  zeros(m)
          zeros(m)     zeros(m) ];
    H = sparse(H);

% ------------------------------------------------------------------
function H = hessian(x, sigma, lambda, auxdata)  
    [n m A y lambda] = deal(auxdata{:});
    H = [ tril(A'*A)  zeros(m)
         zeros(m)   zeros(m) ];
    H = sparse(sigma * H);
  