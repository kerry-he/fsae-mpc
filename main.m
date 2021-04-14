clear all; close all; clc;

%% Add paths to required function folders
addpath(genpath('util'));
addpath(genpath('spline'));
addpath(genpath('mpc'));
addpath(genpath('vehicle_models'));
addpath(genpath('optimizers'));

%% Obtain track spline
filename = "data/fsg2019.csv";
[x, y, vx, vy, ax, ay, dt, rx, ry, lx, ly] = read_raceline_csv(filename);

% Generate spline
x_spline = make_spline_periodic(x);
y_spline = make_spline_periodic(y);
[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 75, true);
kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl); 

%% Set MPC parameters
% Define time horizon
N_x = 5;
N_u = 2;
N_steps = 40;
dt = 0.1;

% Define constraints
state_idx = [4, 5];
soft_idx = [2];
N_soft = length(soft_idx);
x_lb = repmat([0, -0.4, -0.75], N_steps, 1);
x_ub = repmat([1e10, 0.4, 0.75], N_steps, 1);
x_lb = x_lb(:); x_ub = x_ub(:);

u_lb = [repmat([-10; -0.4], N_steps, 1); 0];
u_ub = [repmat([10; 0.4], N_steps, 1); 1e10];

% Define cost weights
Q = [5; 50; 10; 0; 0];
Q_terminal = Q * 10;
R = [10, 10];
R_soft = 1e8;

% Sample parameters
TARGET_VEL = 10;
x_ref = zeros(N_x, N_steps);
x_ref(4, :) = TARGET_VEL;
u_ref = zeros(N_u, N_steps);
x0 = zeros(N_x, 1);

%% Simulate MPC
N_simulation = 500;
x = zeros(4, 1);
x_opt = reshape(x_ref, N_x, N_steps);
x_history = zeros(N_simulation, 4);
u_opt_history = zeros(N_simulation, N_u + N_soft);
x_opt_history = zeros(N_simulation, N_x);

for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    [s, n, mu] = cartesian_to_curvilinear(x(1), x(2), x(3), x_spline, y_spline, dl, x_opt(1));
    x0 = [s; n; mu; x(4); x_opt(5)];
    
    % Define new reference points
    x_ref(1, :) = s : TARGET_VEL*dt : s+TARGET_VEL*dt*(N_steps - 1);
    
    % Define QP problem
    [A, B, d] = linearise_kinematic_curvilinear(reshape(x_opt, N_x, N_steps), u_ref, kappa);
    [A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt);
    [B_bar, xA, lbA, ubA] = state_constraints(A_bar, B_bar, d_bar, x0, x_lb, x_ub, state_idx, soft_idx);
    [H, f] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R, R_soft);
    
    % Solve QP problem
    [u_opt,fval] = qpOASES(H, f, xA, u_lb, u_ub, lbA, ubA);
    x_opt = A_bar*x0 + B_bar*u_opt + d_bar;

    % Update vehicle model
    x = kinematic_bicycle(x, [u_opt(1); x_opt(5)], dt);
    x_history(i, :) = x';
    u_opt_history(i, :) = [u_opt(1:2); u_opt(N_u*N_steps + 1)]';
    x_opt_history(i, :) = x_opt(1:5)';
     
    if mod(i, 50) == 0
        display("Running iteration: " + i)
    end
end

%% Plot results
x_int = interpolate_spline(0:1:L, x_spline, dl);
y_int = interpolate_spline(0:1:L, y_spline, dl);

plot(x_history(:, 1), x_history(:, 2))
hold on
plot(x_int, y_int)
plot(rx, ry, "*")
plot(lx, ly, "*")
