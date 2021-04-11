clear all; close all; clc;

%% Add paths to required function folders
addpath(genpath('util'));
addpath(genpath('spline'));
addpath(genpath('mpc'));
addpath(genpath('vehicle_models'));
addpath(genpath('qpOASES'));

%% Obtain track spline
filename = "data/fsg2019.csv";
[x, y, vx, vy, ax, ay, dt, rx, ry, lx, ly] = read_raceline_csv(filename);

% Generate spline
x_spline = make_spline_periodic(x);
y_spline = make_spline_periodic(y);
[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 75, true);
kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl); 

%% Set MPC parameters
% Define cost weights
N_x = 5;
N_u = 2;
Q = [5; 50; 10; 0; 0];
Q_terminal = Q * 10;
R = [10, 10];
R_soft = [1e8];

% Define time horizon
N_steps = 40;
dt = 0.05;

% Define constraints
u_lb = [repmat([-10; -0.4], N_steps, 1); repmat([-1e10], N_steps, 1)];
u_ub = [repmat([10; 0.4], N_steps, 1); repmat([1e10], N_steps, 1)];

x_lb = repmat([-1.0], N_steps, 1);
x_ub = repmat([1.0], N_steps, 1);

% Sample parameters
TARGET_VEL = 10;
x_ref = zeros(N_x, N_steps);
x_ref(4, :) = TARGET_VEL;
u_ref = zeros(N_u, N_steps);
x0 = zeros(N_x, 1);

%% Simulate MPC
N_simulation = 500;
x = zeros(4, 1);
x_opt = zeros(N_x * N_steps, 1);
u_history = zeros(N_simulation, 2);
x_history = zeros(N_simulation, 4);
u_opt_history = zeros(N_simulation, 3);
x_opt_history = zeros(N_simulation, 5);

for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    s = closest_point(x(1), x(2), x_spline, y_spline, dl, x_opt(1), 0.01);
    
    car_from_track = [x(1) - interpolate_spline(s, x_spline, dl); x(2) - interpolate_spline(s, y_spline, dl)]; 
    tangent = [-interpolate_spline_d(s, y_spline, dl); interpolate_spline_d(s, x_spline, dl)];
    tangent = tangent / norm(tangent);
    n = dot(car_from_track, tangent);
         
    mu = angdiff(interpolate_angle(s, x_spline, y_spline, dl), x(3));
    x0 = [s; n; mu; x(4); x_opt(5)];
    
    % Define new reference points
    x_ref(1, :) = s : TARGET_VEL*dt : s+TARGET_VEL*dt*(N_steps - 1);

    % Define QP problem
    [A, B, d] = linearise_kinematic_curvilinear(x_ref, u_ref, kappa);
    [A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt);
    
    % Soft constraints
    B_bar = [B_bar, zeros(N_x*N_steps, N_steps)];
    [H, f] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R, R_soft);
    
    % Define state constraints
    Al = A_bar(2:N_x:end, :);
    Bl = B_bar(2:N_x:end, :);
    Bl(:, N_u*N_steps + 1 : end) = eye(N_steps);
    Dl = d_bar(2:N_x:end, :);
    offset = Al*x0 + Dl;
    
    lbA = x_lb - offset;
    ubA = x_ub - offset;
    
    % Solve QP problem
    [u_opt,fval] = qpOASES(H, f, Bl, u_lb, u_ub, lbA, ubA);
    x_opt = A_bar*x0 + B_bar*u_opt + d_bar;

    % Update vehicle model
    x = kinematic_bicycle(x, [u_opt(1); x_opt(5)], dt);
    u_history(i, :) = [u_opt(1); x_opt(5)]';
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
