% Auhtor: Rohan Singh

% Testing Gaussian Elimination
B = [3 2 1; 9 8 -6; 3 -5 0; 0 2 6];

%[L1, U1] = luFactorization(B);
%[L2, U2, P2] = luFactorizationPartialPivot(B);

% Matlab LU
[L3, U3] = lu(B);



disp('LU Factorization using MATLAB:');
disp('L:');
disp(L3);
disp('U:');
disp(U3);


% Testing Cholesky
A = [9 -3 3; -3 5 1; 3 1 3];
R = choleskyFactorization(A);

disp('Cholesky Factorization:');
disp('R:');
disp(R);




%% Functions are implemented Here

% Gaussian Elimination without partial pivoting
function [L, U] = luFactorization(A)
    [m, n] = size(A);
    if m ~= n
        error('Input matrix must be square');
    end

    L = eye(n); % Initialize L as identity matrix
    U = A;      % Copy of the input matrix to perform operations

    for k = 1:n-1
        for i = k+1:n
            if U(k,k) == 0
                error('Matrix is singular');
            end
            L(i,k) = U(i,k) / U(k,k);
            U(i,k:n) = U(i,k:n) - L(i,k) * U(k,k:n);
        end
    end
end


% Gaussian Elimination with partial pivoting
function [L, U, P] = luFactorizationPartialPivot(A)
    [m, n] = size(A);
    if m ~= n
        error('Input matrix must be square');
    end

    L = eye(n);
    U = A;
    P = eye(n);

    for k = 1:n-1
        [~, pivot] = max(abs(U(k:n,k))); % Partial pivoting
        pivot = pivot + k - 1;

        if U(pivot, k) == 0
            error('Matrix is singular');
        end

        % Swap rows in U
        temp = U(k,:);
        U(k,:) = U(pivot,:);
        U(pivot,:) = temp;

        % Swap rows in L
        temp = L(k, 1:k-1);
        L(k, 1:k-1) = L(pivot, 1:k-1);
        L(pivot, 1:k-1) = temp;

        %swapping for P
        temp = P(k,:);
        P(k,:) = P(pivot,:);
        P(pivot,:) = temp;

        for i = k+1:n
            L(i,k) = U(i,k) / U(k,k);
            U(i,k:n) = U(i,k:n) - L(i,k) * U(k,k:n);
        end
    end
end


% Function for Cholesky
function R = choleskyFactorization(A)
    [m, n] = size(A);
    if m ~= n || any(any(A ~= A'))
        error('Input matrix must be symmetric');
    end

    L = zeros(n);

    R = A;
    for k = 1:m
        for j = k+1 : m 
            a = R(k,j) / R(k,k);
            R(j,j:m) = R(j,j:m) - R(k,j:m)*a;
        end
        R(k,k:m) = R(k,k:m)/sqrt(R(k,k));
    end
    R = triu(R);

end
