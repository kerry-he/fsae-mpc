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
MODE = "LTV-MPC";

% Define time horizon
N_x = 5;
N_u = 2;
N_steps = 40;
dt = 0.1;

% Sample parameters
TARGET_VEL = 10;
x_ref = zeros(N_x, N_steps);
x_ref(4, :) = TARGET_VEL;
u_ref = zeros(N_u, N_steps);

%% Simulate MPC
N_simulation = 500;
x = zeros(4, 1);
x_opt = reshape(x_ref, N_x, N_steps);
ipopt_info = [];
x0 = zeros(N_x, 1);

x_history = zeros(N_simulation, 4);
u_opt_history = zeros(N_simulation, N_u);
x_opt_history = zeros(N_simulation, N_x);
QP = 0;

for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    [s, n, mu] = cartesian_to_curvilinear(x(1), x(2), x(3), x_spline, y_spline, dl, x_opt(1));
    x0 = [s; n; mu; x(4); x_opt(5)];
    
    % Define new reference points
    x_ref(1, :) = s : TARGET_VEL*dt : s+TARGET_VEL*dt*(N_steps - 1);
    
    % Solve MPC problem
    if MODE == "LTV-MPC"
        % Solve linear time varying MPC problem
        [u_opt, x_opt, QP] = ltvmpc_kinetmatic_curvilinear(x0, x_ref, kappa, dt, ...
            reshape(x_opt, N_x, N_steps), zeros(N_u, N_steps), QP);
    elseif MODE == "NMPC"
        % Solve the nonlinear MPC problem
        [x_mpc, ipopt_info] = nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, ipopt_info);
        x_opt = x_mpc(1:5);
        u_opt = x_mpc(6:7);
    end

    % Update vehicle model
    x = kinematic_bicycle(x, [u_opt(1); x_opt(5)], dt);
    x_history(i, :) = x';
    u_opt_history(i, :) = u_opt(1:2)';
    x_opt_history(i, :) = x_opt(1:5)';
     
    if mod(i, 50) == 0
        display("Running iteration: " + i)
    end
end

if MODE == "LTV-MPC"
	qpOASES_sequence('c', QP);
end

%% Plot results
x_int = interpolate_spline(0:1:L, x_spline, dl);
y_int = interpolate_spline(0:1:L, y_spline, dl);

plot(x_history(:, 1), x_history(:, 2))
hold on
plot(x_int, y_int)
plot(rx, ry, "*")
plot(lx, ly, "*")
