function [x_P_new, y_P_new, dl, L] = arclength_reparam(x_P, y_P, M)
%ARCLENGTH_REPARAM Summary of this function goes here
%   M - number of segments of new spline

    % Calculate total arc length of spline
    N = length(x_P);
    l = zeros(N, 1);
    
    for i = 1:N
        x_d = @(t) -3*(1-t).^2*x_P(i, 1) + 3*(3*t.^2 - 4*t + 1)*x_P(i, 1) ...
            + 3*(2*t - 3*t.^2)*x_P(i, 3) + 3*t.^2*x_P(i, 4);
        y_d = @(t) -3*(1-t).^2*y_P(i, 1) + 3*(3*t.^2 - 4*t + 1)*y_P(i, 1) ...
            + 3*(2*t - 3*t.^2)*y_P(i, 3) + 3*t.^2*y_P(i, 4);
        l(i) = integral(@(t) sqrt(x_d(t).^2 + y_d(t).^2), 0., 1.);
    end
    
    l_cum = cumsum(l);
    l_cum = [0; l_cum];
    dl = l_cum(N + 1) / M;
    
    % Find evenly spaced points along arc
    Px = zeros(M + 1, 1);
    Py = zeros(M + 1, 1);
    Px(1) = x_P(1, 1);
    Py(1) = y_P(1, 1);
    Px(M + 1) = x_P(N, 1);
    Py(M + 1) = y_P(N, 1);
    
    for i = 1:M - 1
        % Define function to solve the root for
        j = find(l_cum >= i*dl, 1) - 1;
        
        x_d = @(t) -3*(1-t).^2*x_P(j, 1) + 3*(3*t.^2 - 4*t + 1)*x_P(j, 1) ...
            + 3*(2*t - 3*t.^2)*x_P(j, 3) + 3*t.^2*x_P(j, 4);
        y_d = @(t) -3*(1-t).^2*y_P(j, 1) + 3*(3*t.^2 - 4*t + 1)*y_P(j, 1) ...
            + 3*(2*t - 3*t.^2)*y_P(j, 3) + 3*t.^2*y_P(j, 4);
        f = @(T) integral(@(t) sqrt(x_d(t).^2 + y_d(t).^2), 0., T) + l_cum(j) - i*dl;
        
        t_i = bisection(0, 1, f, 0.01);
        
        Px(i + 1) = interpolate_spline(t_i + j - 1, x_P, 1.);
        Py(i + 1) = interpolate_spline(t_i + j - 1, y_P, 1.);
    end
    
    % Recalculate spline parameters using new points
    x_P_new = make_spline(Px);
    y_P_new = make_spline(Py);
    
    L = l_cum(N);
    
end

