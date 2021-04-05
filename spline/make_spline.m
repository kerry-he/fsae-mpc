function [coeffs] = make_spline(P)
%MAKE_SPLINE Creates a cubic B-spline with Bezier segments
%   Detailed explanation goes here

    % Define known points
    N = length(P) - 1;
    P0 = P(1:N);
    P3 = P(2:N + 1);

    % Populate matrix A
    diag = zeros(N, 3);
    diag(1, :) = [1, 2, 0];
    diag(2:N - 1, :) = repmat([1, 4, 1], N - 2, 1);
    diag(N - 1, 1) = 2;
    diag(N, :) = [0, 7, 1];
    
    A = spdiags(diag, -1:1, N, N);
    
    % Populate vector b
    b = zeros(N, 1);
    b(1) = P(1) + 2*P(2);
    b(2:N - 1) =  4*P(2:N - 1) + 2*P(3:N);
    b(N) = 8*P(N) + P(N + 1);
    
    % Solve for P1
    P1 = A \ b;
    
    % Backsubstitute to solve for P2
    P2 = zeros(N, 1);
    P2(1) = 2*P1(1) - P(1);
    P2(2:N - 1) = 2*P(3:N) - P1(3:N);
    P2(N) = (P(N + 1) + P1(N)) / 2;
    
    % Output coefficients
    coeffs = [P0, P1, P2, P3];
end

