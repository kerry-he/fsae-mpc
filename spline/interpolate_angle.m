function theta = interpolate_angle(s, x_P, y_P, dl)
%INTERPOLATE_SPLINE Find a value along a spline
%   INPTUS:
%       s - Spline input parameter
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)

%   OUTPUTS:
%       theta - Value of spline at t

    % Precalculate derivatives
    X_d = interpolate_spline_d(s, x_P, dl);
    Y_d = interpolate_spline_d(s, y_P, dl);
    
    % Calculate curvature
    theta = atan2(Y_d, X_d);
end

