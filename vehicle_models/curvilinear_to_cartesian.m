function [x, y, theta] = curvilinear_to_cartesian(s, n, mu, x_P, y_P, dl)
%CARTESIAN_TO_CURVILINEAR Transforms Cartesian position coordinates to
%curvilinear coordinates
%   INPUTS:
%       s - Distance along path
%       n - Distance perpendicular from path
%       mu - Anglear deviation from path tangent
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       x - x Cartesian coordinate
%       y - y Cartesain coordinate
%       theta - Yaw

    x_track = interpolate_spline(s, x_P, dl);
    y_track = interpolate_spline(s, y_P, dl);
    
    x_tangent = -interpolate_spline_d(s, y_P, dl);
    y_tangent = interpolate_spline_d(s, x_P, dl);
    tangent_norm = sqrt(x_tangent.^2 + y_tangent.^2);
    x_tangent = x_tangent ./ tangent_norm;
    y_tangent = y_tangent ./ tangent_norm;
    
    x = x_track + n .* x_tangent;
    y = y_track + n .* y_tangent;
    
    theta = interpolate_angle(s, x_P, y_P, dl) + mu;

end

