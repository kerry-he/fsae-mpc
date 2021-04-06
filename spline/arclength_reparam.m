function [x_P_new, y_P_new, dl, L] = arclength_reparam(x_P, y_P, M, periodic)
%ARCLENGTH_REPARAM Reparametrise a spline to a function of arc length
%   INPUTS:
%       x_P - x spline coefficients
%       y_P - y spline coefficients
%       M - number of segments of new spline
%       periodic - Boolean for if periodic spline should be made
%   OUTPUTS:
%       x_P_new - Reparameterised x spline coefficients
%       y_P_new - Reparameterised y spline coefficients
%       dl - Segment length
%       L - Total length

    % Step 1: Calculate total arc length of spline
    N = length(x_P);
    l = zeros(N, 1);
    
    for i = 1:N
        % Define arc length integral function
        x_d = @(t) -3*(1-t).^2*x_P(i, 1) + 3*(3*t.^2 - 4*t + 1)*x_P(i, 1) ...
            + 3*(2*t - 3*t.^2)*x_P(i, 3) + 3*t.^2*x_P(i, 4);
        y_d = @(t) -3*(1-t).^2*y_P(i, 1) + 3*(3*t.^2 - 4*t + 1)*y_P(i, 1) ...
            + 3*(2*t - 3*t.^2)*y_P(i, 3) + 3*t.^2*y_P(i, 4);
        l(i) = integral(@(t) sqrt(x_d(t).^2 + y_d(t).^2), 0., 1.);
    end
    
    l_cum = cumsum(l);
    l_cum = [0; l_cum];
    dl = l_cum(N + 1) / M;
    
    % Step 2: Find evenly spaced points along arc
    Px = zeros(M + 1, 1);
    Py = zeros(M + 1, 1);
    Px(1) = x_P(1, 1);
    Py(1) = y_P(1, 1);
    Px(M + 1) = x_P(N, 4);
    Py(M + 1) = y_P(N, 4);
    
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
    
    % Step 3: Recalculate spline parameters using new points
    if periodic
        x_P_new = make_spline_periodic(Px(1:M));
        y_P_new = make_spline_periodic(Py(1:M));        
    else
        x_P_new = make_spline(Px);
        y_P_new = make_spline(Py);
    end
    
    L = l_cum(end);
    
end

function x = bisection(xl, xu, f, epsilon)
%BISECTION Bisection method which finds a zero root of a function
%   INPUTS:
%       xl - Initial lower bound
%       xu - Initial upper bound
%       f - Anonymous function to find the root of
%       epsilon - Error threshold
%   OUTPUTS:
%       x - Root of the function

    exit_flag = false;

    while ~exit_flag
        % Bisect between upper and lower bounds
        x = (xl + xu) / 2;
        fx = f(x);
        
        if abs(fx) <= epsilon
            % Within error threshold - exit
            exit_flag = true;
        elseif fx < 0
            % f(x) is negative, so shift lower bound
            xl = x;
        else
            % f(x) is positive, so shift upper bound
            xu = x;
        end
    end
                
end