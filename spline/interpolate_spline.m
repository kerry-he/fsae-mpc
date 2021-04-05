function [x] = interpolate_spline(t, P, dl)
%INTERPOLATE_SPLINE Summary of this function goes here
%   Detailed explanation goes here

    t = t(:);
    i = floor(t / dl) + 1;
    t = t / dl - (i - 1);
    
    x = P(i, 1).*(1-t).^3      + 3*P(i, 2).*(1-t).^2.*t ...
      + 3*P(i, 3).*(1-t).*t.^2 + P(i, 4).*t.^3;

end

