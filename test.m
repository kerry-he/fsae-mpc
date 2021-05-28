clear all; close all; clc;
% Define vehicle constants
m = 200;
I = 200;
lr = 0.6183;
lf = 0.8672;
lr_ratio = lr / (lr + lf);

g = 9.81;

% Define states and controls
s       = 0;
n       = 0;
mu      = 0;
x_d     = 0.01;
y_d     = -5:0.0001:5;
theta_d = 0;
delta   = 0.2; 

% Slip angles
alpha_f = delta - atan((y_d + lf*theta_d) ./ x_d);
alpha_r = -atan((y_d - lr*theta_d) ./ x_d);

% Mass distribution
Fzf = m*g * lf / (lr+lf);
Fzr = m*g * lr / (lr+lf);

% Pacejka magic formula
B = 12.56;
C = 1.38;
D = 1.60;
E = -0.58;

Fcf = Fzf * D * sin(C * atan(B*alpha_f - E*(B*alpha_f - atan(B*alpha_f))));
Fcr = Fzr * D * sin(C * atan(B*alpha_r - E*(B*alpha_r - atan(B*alpha_r))));

Fcf_d = Fzf * D * cos(C * atan(B*alpha_f - E*(B*alpha_f - atan(B*alpha_f)))) ...
            .* C ./ (1 + (B*alpha_f - E*(B*alpha_f - atan(B*alpha_f))).^2) ...
            .* (B - E * (B - B ./ (1 + B^2 * alpha_f.^2)));

Fcr_d = Fzr * D * cos(C * atan(B*alpha_r - E*(B*alpha_r - atan(B*alpha_r)))) ...
            .* C ./ (1 + (B*alpha_r - E*(B*alpha_r - atan(B*alpha_r))).^2) ...
            .* (B - E * (B - B ./ (1 + B^2 * alpha_r.^2)));       
        

% Define common constants
k = 0.1;
denom_nk = 1 ./ (1 - n * k); 
vf = (y_d + lf*theta_d) ./ x_d;
vr = (y_d - lr*theta_d) ./ x_d;
denom_vf2 = 1 ./ (1 + vf.^2);
denom_vr2 = 1 ./ (1 + vr.^2);


Fcr2 = Fcr.^2;
Fcr2d = (Fcr2(2:end) - Fcr2(1:end-1))/0.0001;
Fcr2_d = -2*Fcr.*Fcr_d.*denom_vr2./x_d * (1-exp(-x_d));

% Partial derivatives
s_n = (x_d * cos(mu) - y_d * sin(mu))*denom_nk^2 * k;
s_mu = (-x_d * sin(mu) - y_d * cos(mu))*denom_nk;
s_xd = cos(mu) * denom_nk;
s_yd = -sin(mu) * denom_nk;

n_mu = x_d * cos(mu) - y_d * sin(mu);
n_xd = sin(mu);
n_yd = cos(mu);

mu_n = -s_n * k;
mu_mu = -s_mu * k;
mu_xd = -s_xd * k;
mu_yd = -s_yd * k;
mu_thetad = 1;

xd_xd = -Fcf_d .* denom_vf2 .* vf .* sin(delta) ./ (m * x_d);
xd_yd = (Fcf_d .* denom_vf2 .* sin(delta) ./ x_d + m * theta_d) / m;
xd_thetad = (Fcf_d .* denom_vf2 * lf  .* sin(delta) ./ x_d + m * y_d) / m;
xd_delta = (-Fcf .* cos(delta) - Fcf_d .* sin(delta)) / m;

x = (0 - Fcf.*sin(delta) + m*y_d*theta_d) / m;
xd = (x(2:end) - x(1:end-1))/0.01;

yd_xd = (Fcr_d .* denom_vr2 .* vr ./ x_d + Fcf_d .* denom_vf2 .* vf .* cos(delta) ./ x_d - m * theta_d) / m;
yd_yd = (-Fcr_d  .* denom_vr2 ./ x_d - Fcf_d .* denom_vf2 ./ x_d .* cos(delta)) / m;
yd_thetad = (Fcr_d  .* denom_vr2 .* lr ./ x_d - Fcf_d .* denom_vf2 .* lf ./ x_d .* cos(delta) - m * x_d) / m;
yd_delta = (-Fcf .* sin(delta) + Fcf_d .* cos(delta)) / m;

y = (Fcr + Fcf.*cos(delta) - m*x_d*theta_d) / m;
yd = (y(2:end) - y(1:end-1))/0.001;

t_xd = (lf * Fcf_d .* denom_vf2 .* vf .* cos(delta) ./ x_d - lr * Fcr_d .* denom_vr2 .* vr ./ x_d) / I;
t_yd = (-lf * Fcf_d .* denom_vf2 .* cos(delta) ./ x_d + lr * Fcr_d .* denom_vr2 ./ x_d) / I;
t_thetad = (-lf * Fcf_d .* denom_vf2 * lf .* cos(delta) ./ x_d - lr * Fcr_d .* denom_vr2 .* lr ./ x_d) / I;
t_delta = (-lf * Fcf .* sin(delta) + lf * Fcf_d .* cos(delta)) / I;

t = (lf*Fcf.*cos(delta) - lr*Fcr) / I;
td = (t(2:end) - t(1:end-1))/0.001;
