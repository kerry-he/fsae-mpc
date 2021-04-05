function [x] = bisection(xl, xu, f, epsilon)
%BISECTION Summary of this function goes here
%   Detailed explanation goes here

    exit_flag = false;

    while ~exit_flag
        x = (xl + xu) / 2;
        fx = f(x);
        if abs(fx) <= epsilon
            exit_flag = true;
        elseif fx < 0
            xl = x;
        else
            xu = x;
        end
    end
                
end

