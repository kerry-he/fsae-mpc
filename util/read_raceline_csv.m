function [x, y, vx, vy, ax, ay, dt] = read_raceline_csv(filename)
%READ_RACELINE_CSV Reads an output ORL file and outputs its contents as
%vectors

    % Extract contents into matrix
    A = readmatrix(filename);
    
    % Extract each vector into its own variable
    x = A(:, 1);
    y = A(:, 2);
    vx = A(:, 3);
    vy = A(:, 4);
    ax = A(:, 5);
    ay = A(:, 6);
    dt = A(:, 7);

end

