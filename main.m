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
kappa_d = @(s) interpolate_curvature_d(s, x_spline, y_spline, dl); 

% [x, info, ds, N] = minimum_time_planner(x_spline, y_spline, dl, L);
% [x_pred, y_pred, ~] = curvilinear_to_cartesian(0:ds:ds*(N-1), ...
%     x(1:6:end), x(1:6:end), x_spline, y_spline, dl);

%% Set MPC parameters
MODE = "NMPC";
VISUALISE = true;

% Define time horizon
N_x = 5;
N_u = 2;
N_steps = 20;
dt = 0.1;

% Sample parameters
TARGET_VEL = 20;
x_ref = zeros(N_x, N_steps);
x_ref(4, :) = TARGET_VEL;
u_ref = zeros(N_u, N_steps);

%% Simulate MPC
N_simulation = 350;
x = zeros(5, 1);
x_opt = reshape(x_ref, N_x, N_steps);
u_opt = zeros(N_u*N_steps+1, 1);
x_mpc = repmat([0; 0; 0; 20; 0; 0; 0], N_steps*2+1, 1);
ipopt_info = [];
x0 = zeros(N_x, 1);

x_history = zeros(N_simulation, 5);
u_opt_history = zeros(N_simulation, N_u);
x_opt_history = zeros(N_simulation, N_x);
QP = 0;

figure(1)
plot(rx, ry, "y*")
hold on
plot(lx, ly, "b*")


for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    [s, n, mu] = cartesian_to_curvilinear(x(1), x(2), x(3), x_spline, y_spline, dl, x_opt(1));
    x0 = [s; n; mu; x(4); x(5)];
        
    % Define new reference points
    x_ref(1, :) = x0(1)+TARGET_VEL*dt : TARGET_VEL*dt : x0(1)+TARGET_VEL*dt*N_steps;
    
    % Solve MPC problem
    if MODE == "LTV-MPC"
        % Solve linear time varying MPC problem
        [u_opt, x_opt, QP] = ltvmpc_kinetmatic_curvilinear(x0, x_ref, kappa, kappa_d, dt, ...
            reshape(x_opt, N_x, N_steps), reshape(u_opt(1:end-1), N_u, N_steps), QP);
    elseif MODE == "NMPC"
        % Solve the nonlinear MPC problem
        [x_mpc, ipopt_info] = hs_nmpc_kinematic_curvilinear(x0, x_ref, kappa, kappa_d, dt, x_mpc, ipopt_info);
        x_opt = x_mpc([1:7:end; 2:7:end; 3:7:end; 4:7:end; 5:7:end]);
        x_opt = x_opt(:);
        u_opt = x_mpc([6:7:end; 7:7:end;]);
        u_opt = u_opt(:);
    end

    if VISUALISE
        visualise_mpc(x, x_opt, u_opt, x_spline, y_spline, dl, dt/2)
    end

    [x_pred, y_pred, ~] = curvilinear_to_cartesian(x_opt(1:N_x:end), ...
        x_opt(2:N_x:end), x_opt(3:N_x:end), x_spline, y_spline, dl);
    x_cart_pred = kinematic_bicycle_horizon(x, [u_opt(1:2:N_u*N_steps), ...
        u_opt(2:2:N_u*N_steps)]', dt);
    
    cpu_time(i) = ipopt_info.cpu;
%     error(i) = norm([x_pred, y_pred]' - x_cart_pred(1:2, 2:end));
    
    % Update vehicle model
    x = kinematic_bicycle(x, [u_opt(3); u_opt(4)], dt/2);
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
figure(1)
x_int = interpolate_spline(0:1:L, x_spline, dl);
y_int = interpolate_spline(0:1:L, y_spline, dl);

plot(x_history(:, 1), x_history(:, 2), 'r--')
hold on
plot(x_int, y_int, 'k')
plot(rx, ry, "y*")
plot(lx, ly, "b*")
