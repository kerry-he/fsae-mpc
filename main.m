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

% [x_opt_traj, t, info, ds, N_s, slack] = dynamic_minimum_time_planner(x_spline, y_spline, dl, L);
% [x_pred, y_pred, ~] = curvilinear_to_cartesian(0:ds:ds*(N_s-1), ...
%     x_opt_traj(1:8:end), x_opt_traj(2:8:end), x_spline, y_spline, dl);

%% Setup MPC parameters
MODE = "LTV-MPC";
MODEL = "DYNAMIC";
VISUALISE = true;

% Define time horizon
if MODEL == "KINEMATIC"
    N_x = 5;
elseif MODEL == "DYNAMIC"
    N_x = 7;
end
N_u = 2;
N_steps = 40;
dt = 0.05;

% Initialise reference trajectory as constant 
TARGET_VEL = 20;
x_ref = zeros(N_x, N_steps);
u_ref = zeros(N_u, N_steps);

% Define MPC initial guess
x_mpc = [zeros(N_x, N_steps); zeros(N_u, N_steps)];
x_mpc(1, :) = 10 * (dt:dt:dt*N_steps).^2/2; % Quadratic arc length 
x_mpc(4, :) = 10*dt:10*dt:10*dt*N_steps; % Linear velocity
x_mpc(end-1, :) = 10; % Constant acceleration

if MODE == "C-NMPC"
    x_mpc = [x_mpc, x_mpc(:, end)];
end

x_opt = x_mpc(1:N_x, :);
u_opt = x_mpc(end-N_u+1:end, :);

ipopt_info = [];
QP = 0;

%% Setup simulation parameters
% Define simulation vehicle model initial conditions
N_simulation = 1000;
x = zeros(7, 1);

% Preallocate metric vectors
n_list = zeros(N_simulation, 1);
exit_status = zeros(N_simulation, 1);
cpu_time = zeros(N_simulation, 1);
objective = zeros(N_simulation, 1);
slack_n = zeros(N_simulation, 1);
slack_tyre = zeros(N_simulation, 1);

x_history = zeros(N_simulation, 7);
u_opt_history = zeros(N_simulation, N_u);
x_opt_history = zeros(N_simulation, N_x);

% Initialise map plot
figure(1)
plot(rx, ry, "y*")
hold on
plot(lx, ly, "b*")

%% Setup actuator controllers
vel_pid_settings = {16000.0, 0, 0, 2800};
vel_pid_status = {0, 0};

steer_pid_settings = {80.0, 0, 0, 0.8};
steer_pid_status = {0, 0};

%% Simulate MPC
for i = 1:N_simulation
    % Calculate coordinates in curvilinear frame
    [s, n, mu] = cartesian_to_curvilinear(x(1), x(2), x(3), x_spline, y_spline, dl, x_opt(1));
    if MODEL == "KINEMATIC"
        x0 = [s; n; mu; norm(x(4:5)) ; x(7)];
    elseif MODEL == "DYNAMIC"
        x0 = [s; n; mu; x(4); x(5); x(6); x(7)];
    end
    
    % If finished lap, terminate simulation
    n_list(i) = n;
    if s >= L
        break
    end
    
    % Define new reference points
    if x(4) < TARGET_VEL
        x_ref(4, :) = x0(4)+10*dt : 10*dt : x0(4)+10*dt*N_steps;
        x_ref(4, :) = min(x_ref(4, :), TARGET_VEL);
    else
        x_ref(4, :) = x0(4)-10*dt : -10*dt : x0(4)-10*dt*N_steps;
        x_ref(4, :) = max(x_ref(4, :), TARGET_VEL);
    end
    x_ref(1, :) = x0(1) + cumsum(x_ref(4, :)*dt);
%     x_ref = obtain_reference(x_opt_traj, ds, N_s, t, s, dt, N_steps);
    
    % Solve MPC problem
    if MODE == "LTV-MPC"
        % Solve linear time varying MPC problem
        tic
        if MODEL == "KINEMATIC"
            [u_opt, x_opt, QP, exitflag, fval, slack_opt] = ltvmpc_kinetmatic_curvilinear(x0, x_ref, kappa, dt, ...
                reshape(x_opt, N_x, N_steps), reshape(u_opt, N_u, N_steps), QP);
        elseif MODEL == "DYNAMIC"
            [u_opt, x_opt, QP, exitflag, fval, slack_opt] = ltvmpc_dynamic_curvilinear(x0, x_ref, kappa, dt, ...
            	reshape(x_opt, N_x, N_steps), reshape(u_opt, N_u, N_steps), QP);
        end
        
        exit_status(i) = exitflag;
        objective(i) = fval;        
        cpu_time(i) = toc;
        slack_n(i) = slack_opt(1);
        slack_tyre(i) = slack_opt(4);
        
    else
        % Solve the nonlinear MPC problem
        if MODEL == "KINEMATIC"
            if MODE == "MS-NMPC"
                [x_mpc, slack_mpc, ipopt_info] = rk2_nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, x_mpc(:), ipopt_info);
            elseif MODE == "C-NMPC"
                [x_mpc, slack_mpc, ipopt_info] = trapezoidal_nmpc_kinematic_curvilinear(x0, x_ref, kappa, dt, x_mpc(:), ipopt_info);
            end
            x_opt = x_mpc([1:7:end; 2:7:end; 3:7:end; 4:7:end; 5:7:end]); x_opt = x_opt(:);
            u_opt = x_mpc([6:7:end; 7:7:end;]); u_opt = u_opt(:);            
        elseif MODEL == "DYNAMIC"
            if MODE == "MS-NMPC"
                [x_mpc, slack_mpc, ipopt_info] = rk2_nmpc_dynamic_curvilinear(x0, x_ref, kappa, dt, x_mpc(:), ipopt_info);
            elseif MODE == "C-NMPC"
                [x_mpc, slack_mpc, ipopt_info] = trapezoidal_nmpc_dynamic_curvilinear(x0, x_ref, kappa, dt, x_mpc(:), ipopt_info);
            end            
            x_opt = x_mpc([1:9:end; 2:9:end; 3:9:end; 4:9:end; 5:9:end; 6:9:end; 7:9:end]); x_opt = x_opt(:);
            u_opt = x_mpc([8:9:end; 9:9:end;]); u_opt = u_opt(:);
        end
        
        exit_status(i) = ipopt_info.status;
        objective(i) = ipopt_info.objective;
        cpu_time(i) = ipopt_info.cpu;
        slack_n(i) = slack_mpc(1);
        slack_tyre(i) = slack_mpc(1);
    end

    if VISUALISE
        visualise_mpc(x, x_opt, x_spline, y_spline, dl, MODEL)
    end
    
    % Update vehicle model
    if MODE == "C-NMPC"
        v_ref = x_opt(4 + N_x);
        delta_ref = x_opt(2 * N_x);
    else
        v_ref = x_opt(4);
        delta_ref = x_opt(N_x);
    end
    
    for j = 1:10
        [vel_rate, vel_pid_status] = pid_controller(v_ref, x(4), vel_pid_settings, vel_pid_status);
        [steer_rate, steer_pid_status] = pid_controller(delta_ref, x(7), steer_pid_settings, steer_pid_status);
        x = integrate_cart_dyn(x, [vel_rate; steer_rate], dt/10);
    end
    [~, Fcr] = f_curv_dyn(x, [vel_rate; steer_rate], kappa);
    Fcr_list(i) = Fcr;
    Fx_list(i) = u_opt(1);
    x_history(i, :) = x';
    u_opt_history(i, :) = u_opt(1:2)';
    x_opt_history(i, :) = x_opt(1:N_x)';

    if mod(i, 50) == 0
        display("Running iteration: " + i)
    end
end

if MODE == "LTV-MPC"
	qpOASES_sequence('c', QP);
end

%% Metrics

slack_all = (slack_n(1:i-1)==0) & (slack_tyre(1:i-1)==0);
friction_ellipse = (Fcr_list / (280*9.163)).^2 + (Fx_list / 10.0).^2;

COPY_FORMAT = true; 
if COPY_FORMAT
    fprintf('%f\n', (i-1)*dt)
    fprintf('%f\n', sum(abs(n_list(abs(n_list)>0.75)) - 0.75) * dt)
    fprintf('%f\n', max(abs(n_list(abs(n_list)>0.75)) - 0.75))
    fprintf('%.8f\n', mean(cpu_time(1:i-1)))
    fprintf('%.8f\n', median(cpu_time(1:i-1)))
    fprintf('%.8f\n', max(cpu_time(1:i-1)))
    fprintf('%.8f%%\n', sum(exit_status(1:i-1)~=0)/(i-1)*100)
    fprintf('%.8f\n', mean(objective(slack_all)))
    fprintf('%.8f%%\n', sum(slack_n(1:i-1)~=0)/(i-1)*100)
    fprintf('%f\n', sum(friction_ellipse(friction_ellipse>1.0) - 1.0) * dt)
    fprintf('%f\n', max(friction_ellipse(friction_ellipse>1.0) - 1.0))
    fprintf('%.8f%%\n', sum(slack_tyre(1:i-1)~=0)/(i-1)*100)
else
    fprintf('Lap time: %f\n', (i-1)*dt)
    fprintf('Track violation: %f\n', sum(abs(n_list(abs(n_list)>0.75)) - 0.75) * dt)
    fprintf('Max track violation: %f\n', max(abs(n_list(abs(n_list)>0.75)) - 0.75))
    fprintf('Average CPU time: %.8f\n', mean(cpu_time(1:i-1)))
    fprintf('Mode CPU time: %.8f\n', median(cpu_time(1:i-1)))
    fprintf('Max CPU time: %.8f\n', max(cpu_time(1:i-1)))
    fprintf('Abnormal exits: %.8f%%\n', sum(exit_status(1:i-1)~=0)/(i-1)*100)
    fprintf('Optimal value: %.8f\n', mean(objective(slack_all)))
    fprintf('Track slack violations: %.8f%%\n', sum(slack_n(1:i-1)~=0)/(i-1)*100)
    fprintf('Tyre violation: %f\n', sum(friction_ellipse(friction_ellipse>1.0) - 1.0) * dt)
    fprintf('Max tyre violation: %f\n', max(friction_ellipse(friction_ellipse>1.0) - 1.0))
    fprintf('Tyre slack violations: %.8f%%\n', sum(slack_tyre(1:i-1)~=0)/(i-1)*100)
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
