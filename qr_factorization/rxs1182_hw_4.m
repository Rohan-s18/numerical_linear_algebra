% Author: Rohan Singh
% Matlab code for HW4 for MATH 431
% This script only contains code for Question 2 of the homework, the house-
% holder code is in the other script

%% (a) Simulating Matrix A as described in the book
[U,X] = qr(randn(80));
[V,X] = qr(randn(80));
S = diag(2.^(-1:-1:-80));

A = U*S*V;

%% (b) Running both QR Factorizations
[QC, RC] = clgs(A);
[QM, RM] = mgs(A);

%% (c) Plotting the values of j and r_jj
rjj_classical = diag(RC);

% Perform modified Gram-Schmidt orthogonalization
rjj_modified = diag(RM);

% Plot the rjj values for classical and modified Gram-Schmidt
j_values = 1:80;
%ylim([1e-25, 1e0]);
xlim([1, 80]);
figure;
scatter(j_values, rjj_classical, 80, 'b','filled', 'DisplayName', 'Classical Gram-Schmidt');
hold on;
scatter(j_values, rjj_modified, 80, 'r', 'filled','DisplayName', 'Modified Gram-Schmidt');
xlabel('j');
ylabel('r_{jj}');
title('r_{jj} Values for Classical vs. Modified Gram-Schmidt');
legend('Location', 'Best');
grid on;


%% Function for Classical Gram Schmitt Orthogonalization
function [Q,R] = clgs(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end

%% Function for Modified Gram Schmitt Orthogonalization
function [Q,R] = mgs(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);

    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * v;
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end
































