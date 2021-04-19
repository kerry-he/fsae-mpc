function x_ddd = interpolate_spline_ddd(t, P, dl)
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
    
	% Spline function derivative equation
    x_ddd = -6*P(i, 1) + 18*P(i, 2) - 18*P(i, 3) + 6*P(i, 4);
   
	% Chain rule
	x_ddd = x_ddd / dl^3;

end

