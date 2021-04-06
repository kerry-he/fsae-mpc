function coeffs = make_spline_periodic(P)
%MAKE_SPLINE Creates a periodic cubic B-spline with Bezier segments
%   INPUTS:
%       P - Input points to fit a curve to (1D)
%   OUTPUTS:
%       coeffs - Coefficients of the spline [P0, P1, P2, P3]

    % Define known points
    N = length(P);
    P0 = P(1:N);
    P3 = P([2:N, 1]);

    % Populate matrix A
    diag = repmat([1, 4, 1], N, 1);
    A = spdiags(diag, -1:1, N, N);
    A(N, 1) = 1;
    A(1, N) = 1;

    % Populate vector b
    b = zeros(N, 1);
    b(1:N - 1) =  4*P(1:N - 1) + 2*P(2:N);
    b(N) = 4*P(N) + 2*P(1);
    
    % Solve for P1
    P1 = A \ b;
    
    % Backsubstitute to solve for P2
    P2 = zeros(N, 1);
    P2(1:N - 1) = 2*P(2:N) - P1(2:N);
    P2(N) = 2*P(1) - P1(1);
    
    % Output coefficients
    coeffs = [P0, P1, P2, P3];
end

