clear all; close all; clc;

%% Add paths to required function folders
addpath(genpath('util'));
addpath(genpath('spline'));

%% Read data
filename = "data/fsg2019.csv";
[x, y, vx, vy, ax, ay, dt] = read_raceline_csv(filename);

x_spline = make_spline_periodic(x);
y_spline = make_spline_periodic(y);
[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 40, true);

x_int = interpolate_spline(0:1:L, x_spline, dl);
y_int = interpolate_spline(0:1:L, y_spline, dl);

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