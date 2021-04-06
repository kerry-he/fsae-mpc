function kappa = interpolate_curvature(s, x_P, y_P, dl)
%INTERPOLATE_CURVATURE Find the curvature along a spline
%   INPTUS:
%       s - Spline input parameter
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       kappa - Curvature of spline at point

    % Precalculate derivatives
    X_d = interpolate_spline_d(s, x_P, dl);
    Y_d = interpolate_spline_d(s, y_P, dl);
    X_dd = interpolate_spline_dd(s, x_P, dl);
    Y_dd = interpolate_spline_dd(s, y_P, dl);
    
    % Calculate curvature
    kappa = (X_d.*Y_dd - X_dd.*Y_d) ./ (X_d.^2 + Y_d.^2).^(3/2);

end

