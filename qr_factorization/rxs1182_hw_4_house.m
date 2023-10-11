% Author: Rohan Singh
% Matlab code for HW4 for MATH 431
% This script only contains code for Question 3 of the homework (the house-
% holder code), the code for Gram-Schmitt is in the other script

%% Example Run for 10.2
A = rand(4, 3);
[W, R] = house(A);
Q = formQ(W);

disp('Question 10.2');
disp('A:');
disp(A);
disp('QR:');
disp(Q*R);


%% Question 10.3 Running the example for QR factorization of Z
Z = [1, 2, 3;
     4, 5, 6;
     7, 8, 7;
     4, 2, 3;
     4, 2, 2];

% Running modified gram schmitt
[Q_mgs, R_mgs] = mgs(Z);

% Runnning the householder routines
[W_house, R_house] = house(Z);
Q_house = formQ(W_house);

% Running the matlab version
[Q_mat, R_mat] = qr(Z,0);

disp('Question 10.3');
disp('Modified Gram Schmitt Q and R:');
disp(Q_mgs);
disp(R_mgs);

disp('\n Householder Q and R:');
disp(Q_house);
disp(R_house);

disp('\n Matlab Q and R:');
disp(Q_mat);
disp(R_mat);





%% Question 10.2 (a) Code for implicit QR factorization using householder
function [W, R] = house(A)
    [m, n] = size(A);
    W = zeros(m, n);
    R = A;

    for k = 1:n
        x = R(k:m, k);
        v = zeros(m - k + 1, 1);
        v(1) = sign(x(1)) * norm(x) + x(1);
        v = v / norm(v);
        R(k:m, k:n) = R(k:m, k:n) - 2 * v * (v' * R(k:m, k:n));
        W(k:m, k) = v;
    end
end

%% Question 10.2 (b) Code to form Q from W
function Q = formQ(W)
    [m, n] = size(W);
    Q = eye(m);

    for k = n:-1:1
        v = W(k:m, k);
        Q(k:m, :) = Q(k:m, :) - 2 * v * (v' * Q(k:m, :));
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




