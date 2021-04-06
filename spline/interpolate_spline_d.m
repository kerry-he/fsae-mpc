function x_d = interpolate_spline_d(t, P, dl)
%INTERPOLATE_SPLINE_D Find the derivative along a spline
%   INPTUS:
%       t - Spline input parameter
%       P - Spline coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       x_d - Derivative of spline at t

    % Find segment corresponding to input
    t = t(:); % Ensure t is a column vector
    t = mod(t, dl*length(P));
    i = floor(t / dl) + 1;
    t = t / dl - (i - 1);
    
    % Spline function derivative equation
    x_d = -3*(1-t).^2.*P(i, 1)       + 3*(3*t.^2 - 4*t + 1).*P(i, 2) ...
       + 3*(2*t - 3*t.^2).*P(i, 3) + 3*t.^2.*P(i, 4);
   
	% Chain rule
	x_d = x_d / dl;

end

