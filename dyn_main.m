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
[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 100, true);
kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl); 

% [x, info, ds, N] = minimum_time_planner(x_spline, y_spline, dl, L);
% [x_pred, y_pred, ~] = curvilinear_to_cartesian(0:ds:ds*(N-1), ...
%     x(1:6:end), x(1:6:end), x_spline, y_spline, dl);

%% Set MPC parameters
MODE = "NMPC";
VISUALISE = true;

% Define time horizon
N_x = 7;
N_u = 2;
N_steps = 40;
dt = 0.05;

% Sample parameters
TARGET_VEL = 20;
x_ref = zeros(N_x, N_steps);
x_ref(4, :) = TARGET_VEL;
u_ref = zeros(N_u, N_steps);

%% Set up actuator controllers
vel_pid_settings = {16000.0, 0, 0, 2000};
vel_pid_status = {0, 0};

steer_pid_settings = {80.0, 0, 0, 0.8};
steer_pid_status = {0, 0};

%% Simulate MPC
N_simulation = 1000;
x = zeros(7, 1);
x(4) = 10;
x_opt = reshape(x_ref, N_x, N_steps);
u_opt = zeros(N_u*N_steps, 1);
x_mpc = [x_opt; zeros(N_u, N_steps)];
% x_mpc = repmat([0; 0; 0; 20; 0; 0; 0], N_steps*2 + 1, 1);
x_mpc = [x_mpc(:); 0];
ipopt_info = [];
x0 = zeros(N_x, 1);

n_list = zeros(N_simulation, 1);
exit_status = zeros(N_simulation, 1);
cpu_time = zeros(N_simulation, 1);
objective = zeros(N_simulation, 1);
slack = zeros(N_simulation, 1);

x_history = zeros(N_simulation, 7);
u_opt_history = zeros(N_simulation, N_u);
x_opt_history = zeros(N_simulation, N_x);
QP = 0;

figure(1)
plot(rx, ry, "y*")
hold on
plot(lx, ly, "b*")
temp=[];

for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    [s, n, mu] = cartesian_to_curvilinear(x(1), x(2), x(3), x_spline, y_spline, dl, x_opt(1));
    x0 = [s; n; mu; x(4); x(5); x(6); x(7)];
    
    n_list(i) = n;
    if s >= L
        break
    end
        
    % Define new reference points
    x_ref(1, :) = x0(1)+TARGET_VEL*dt : TARGET_VEL*dt : x0(1)+TARGET_VEL*dt*N_steps;
    
    % Solve MPC problem
    if MODE == "LTV-MPC"
        % Solve linear time varying MPC problem
        tic
        [u_opt, x_opt, QP, exitflag, fval] = ltvmpc_kinetmatic_curvilinear(x0, x_ref, kappa, dt, ...
            reshape(x_opt, N_x, N_steps), reshape(u_opt(1:end-1), N_u, N_steps), QP);
        
        exit_status(i) = exitflag;
        objective(i) = fval;        
        cpu_time(i) = toc;
        slack(i) = u_opt(end);
        
    elseif MODE == "NMPC"
        % Solve the nonlinear MPC problem
        [x_mpc, ipopt_info] = euler_nmpc_dynamic_curvilinear(x0, x_ref, kappa, dt, x_mpc, ipopt_info);
        x_opt = x_mpc([1:9:end-1; 2:9:end-1; 3:9:end-1; 4:9:end-1; 5:9:end-1; 6:9:end-1; 7:9:end-1]);
        x_opt = x_opt(:);
        u_opt = x_mpc([8:9:end-1; 9:9:end-1;]);
        u_opt = u_opt(:);
        
        exit_status(i) = ipopt_info.status;
        objective(i) = ipopt_info.objective;
        cpu_time(i) = ipopt_info.cpu;
        slack(i) = x_mpc(end);
    end

    if VISUALISE
        visualise_mpc(x, x_opt, u_opt, x_spline, y_spline, dl, dt)
    end
    
    % Update vehicle model
    for j = 1:10
        [vel_rate, vel_pid_status] = pid_controller(x_opt(4), x(4), vel_pid_settings, vel_pid_status);
        [steer_rate, steer_pid_status] = pid_controller(x_opt(7), x(7), steer_pid_settings, steer_pid_status);
        x = integrate_cart_dyn(x, [vel_rate; steer_rate], dt/10);
    end
    x_history(i, :) = x';
    u_opt_history(i, :) = u_opt(1:2)';
    x_opt_history(i, :) = x_opt(1:7)';
    
    if mod(i, 50) == 0
        display("Running iteration: " + i)
    end
end

if MODE == "LTV-MPC"
	qpOASES_sequence('c', QP);
end

%% Metrics
COPY_FORMAT = false; 
if COPY_FORMAT
    fprintf('%f\n', (i-1)*dt)
    fprintf('%f\n', sum(abs(n_list(abs(n_list)>0.75)) - 0.75) * dt)
    fprintf('%f\n', max(abs(n_list(abs(n_list)>0.75)) - 0.75))
    fprintf('%.8f\n', mean(cpu_time(1:i-1)))
    fprintf('%.8f\n', max(cpu_time(1:i-1)))
    fprintf('%.8f%%\n', sum(exit_status(1:i-1)~=0)/(i-1)*100)
    fprintf('%.8f\n', mean(objective(slack(1:i-1)==0)))
    fprintf('%.8f%%\n', sum(slack(1:i-1)~=0)/(i-1)*100)
else
    fprintf('Lap time: %f\n', (i-1)*dt)
    fprintf('Track violation: %f\n', sum(abs(n_list(abs(n_list)>0.75)) - 0.75) * dt)
    fprintf('Max track violation: %f\n', max(abs(n_list(abs(n_list)>0.75)) - 0.75))
    fprintf('Average CPU time: %.8f\n', mean(cpu_time(1:i-1)))
    fprintf('Max CPU time: %.8f\n', max(cpu_time(1:i-1)))
    fprintf('Abnormal exits: %.8f%%\n', sum(exit_status(1:i-1)~=0)/(i-1)*100)
    fprintf('Optimal value: %.8f\n', mean(objective(slack(1:i-1)==0)))
    fprintf('Slack violations: %.8f%%\n', sum(slack(1:i-1)~=0)/(i-1)*100)    
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
