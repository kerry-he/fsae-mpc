function x = interpolate_spline(t, P, dl)
%INTERPOLATE_SPLINE Find a value along a spline
%   INPTUS:
%       t - Spline input parameter
%       P - Spline coefficients
%       dl - Scale factor (for arclength parameterised spline, otherwise =1)
%   OUTPUTS:
%       x - Value of spline at t

    % Find segment corresponding to input
    t = t(:); % Ensure t is a column vector
    t = mod(t, dl*length(P));
    i = floor(t / dl) + 1;
    t = t / dl - (i - 1);
    
    % Spline function equation
    x = P(i, 1).*(1-t).^3      + 3*P(i, 2).*(1-t).^2.*t ...
      + 3*P(i, 3).*(1-t).*t.^2 + P(i, 4).*t.^3;

end

