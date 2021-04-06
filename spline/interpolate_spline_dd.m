function x_dd = interpolate_spline_dd(t, P, dl)
%INTERPOLATE_SPLINE_DD Find the second derivative along a spline
%   INPTUS:
%       t - Spline input parameter
%       P - Spline coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       x_dd - Second derivative of spline at t

    % Find segment corresponding to input
    t = t(:); % Ensure t is a column vector
    t = mod(t, dl*length(P));
    i = floor(t / dl) + 1;
    t = t / dl - (i - 1);
    
	% Spline function derivative equation
    x_dd = 6*(1-t).*P(i, 1)     + 6*(3*t - 2).*P(i, 2) ...
         + 6*(1 - 3*t).*P(i, 3) + 6*t.*P(i, 4);
   
	% Chain rule
	x_dd = x_dd / dl^2;

end

