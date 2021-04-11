function [s, n, mu] = cartesian_to_curvilinear(x, y, theta, x_P, y_P, dl, s0)
%CARTESIAN_TO_CURVILINEAR Transforms Cartesian position coordinates to
%curvilinear coordinates
%   INPUTS:
%       x - x Cartesian coordinate
%       y - y Cartesain coordinate
%       theta - Yaw
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%       s0 - Initial guess for position
%   OUTPUTS:
%       s - Distance along path
%       n - Distance perpendicular from path
%       mu - Anglear deviation from path tangent

    s = closest_point(x, y, x_P, y_P, dl, s0, 0.01);
    
    car_from_track = [x - interpolate_spline(s, x_P, dl); ...
                      y - interpolate_spline(s, y_P, dl)]; 
    tangent = [-interpolate_spline_d(s, y_P, dl); ...
                interpolate_spline_d(s, x_P, dl)];
    tangent = tangent / norm(tangent);
    n = dot(car_from_track, tangent);

    mu = angdiff(interpolate_angle(s, x_P, y_P, dl), theta);

end

