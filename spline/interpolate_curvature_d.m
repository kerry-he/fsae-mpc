function kappa_d = interpolate_curvature_d(s, x_P, y_P, dl)
%INTERPOLATE_CURVATURE Find the curvature along a spline
%   INPTUS:
%       s - Spline input parameter
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       kappa - Curvature of spline at point

    % Finite differenes
    delta = dl/10;
    kappa_l = interpolate_curvature(s - delta, x_P, y_P, dl);
    kappa_u = interpolate_curvature(s + delta, x_P, y_P, dl);
    
    % Calculate curvature
    kappa_d = (kappa_u - kappa_l) / (2 * delta);

end

