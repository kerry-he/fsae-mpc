function s = closest_point(x0, y0, x_P, y_P, dl, s, epsilon)
%CLOSEST_POINT_SPLINE Find the point on a spline that is closest to a given
%point using Newton-Raphson root finding optimization
%   INPUT: 
%       x0 - Input point x coordinate
%       y0 - Input point y cooridnate
%       x_P - Spline x coefficients
%       y_P - Spline y coefficients
%       dl - Segment length
%       s - Initial guess
%       epsilon - Error threshold
%   OUTPUTS:
%       s - Parameter representing closest point to spline

    delta = epsilon * 2;

    while abs(delta) > epsilon
        % Calculate derivatives required for Newton step
        X = interpolate_spline(s, x_P, dl);
        Y = interpolate_spline(s, y_P, dl);
        X_d = interpolate_spline_d(s, x_P, dl);
        Y_d = interpolate_spline_d(s, y_P, dl);
        X_dd = interpolate_spline_dd(s, x_P, dl);
        Y_dd = interpolate_spline_dd(s, y_P, dl);
        
        dist_d = 2*(X-x0)*X_d + 2*(Y-y0)*Y_d;
        dist_dd = 2*(X-x0)*X_dd + 2*X_d^2 + 2*(Y-y0)*Y_dd + 2*Y_d^2;
        
        % Perform Newton-Raphson step
        delta = dist_d / dist_dd;
        s = s - delta;
    end

end

