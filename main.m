clear all; close all; clc;

%% Add paths to required function folders
addpath(genpath('util'));
addpath(genpath('spline'));
addpath(genpath('mpc'));

%% Obtain track spline
filename = "data/fsg2019.csv";
[x, y, vx, vy, ax, ay, dt] = read_raceline_csv(filename);

% Generate spline
x_spline = make_spline_periodic(x);
y_spline = make_spline_periodic(y);
[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 40, true);
kappa = @(s) interpolate_curvature(s, x_spline, y_spline, dl); 

%% Set MPC parameters
% Define cost weights
Q = [10; 10; 10; 10; 0];
Q_terminal = [50; 50; 50; 50; 0];
R = [10, 10];

% Define time horizon
N_steps = 20;
dt = 0.01;

% Sample parameters
x_ref = rand(5, N_steps);
u_ref = zeros(2, N_steps);
x0 = zeros(5, 1);

% Solve for QP
[A, B, d] = linearise_kinematic_curvilinear(x_ref, u_ref, kappa);
[A_bar, B_bar, d_bar] = sequential_integration(A, B, d, dt);
[H, f] = generate_qp(A_bar, B_bar, d_bar, x0, x_ref, Q, Q_terminal, R);

%% Junk
% x_int = interpolate_spline(0:1:L, x_spline, dl);
% y_int = interpolate_spline(0:1:L, y_spline, dl);
%
% s = closest_point(-5, -5, x_spline, y_spline, dl, 0, 0.01);
% 
% x_c = interpolate_spline(s, x_spline, dl);
% y_c = interpolate_spline(s, y_spline, dl);
% 
% kappa = interpolate_curvature(0:1:L, x_spline, y_spline, dl);
% store = zeros(101, 25);
% for i = 1:25
%     % Define function to solve the root for
%     x_d = @(t) -3*(1-t).^2*x_spline(i, 1) + 3*(3*t.^2 - 4*t + 1)*x_spline(i, 1) ...
%         + 3*(2*t - 3*t.^2)*x_spline(i, 3) + 3*t.^2*x_spline(i, 4);
%     y_d = @(t) -3*(1-t).^2*y_spline(i, 1) + 3*(3*t.^2 - 4*t + 1)*y_spline(i, 1) ...
%         + 3*(2*t - 3*t.^2)*y_spline(i, 3) + 3*t.^2*y_spline(i, 4);
%     f = @(T) integral(@(t) sqrt(x_d(t).^2 + y_d(t).^2), 0., T);
% 
%     in = 0:0.01:1;
%     for j = 1:101
%         store(j, i) = f(in(j));
%     end
% end