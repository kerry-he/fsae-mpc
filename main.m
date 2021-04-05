clear all; close all; clc;

%% Add paths to required function folders
addpath(genpath('util'));
addpath(genpath('spline'));

%% Read data
filename = "data/fsg2019.csv";
[x, y, vx, vy, ax, ay, dt] = read_raceline_csv(filename);

x_spline = make_spline(x);
y_spline = make_spline(y);

x_int = interpolate_spline(0:0.1:75.9, x_spline, 1.0);
y_int = interpolate_spline(0:0.1:75.9, y_spline, 1.0);

[x_spline, y_spline, dl, L] = arclength_reparam(x_spline, y_spline, 50);

x_int = interpolate_spline(0:10:L-0.1, x_spline, dl);
y_int = interpolate_spline(0:10:L-0.1, y_spline, dl);