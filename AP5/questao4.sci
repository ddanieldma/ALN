function [U, R] = qr_House(A)
    // Matriz que a calcula a decomposição QR de A usando os
    // Refletores de Householder.

    [m, n] = size(A)
    
    k = min(m, n)
    R = zeros(m, n)
    U = zeros(m, n)

    for j = 1:k
        x = A(j:m, j)
        x_norm = norm(x, 2)
        if x(1) < 0 then
            x(1) = x(1) - x_norm
        else
            x(1) = x(1) + x_norm
        end

        u = x/x_norm
        U(j:m, j) = u
        A(j:m, j:n) = A(j:m, j:n) - 2 * u * (u' * A(j:m, j:n))
    end

    R = triu(A)
endfunction

function [U, R] = qr_House(A)
    // Matriz que a calcula a decomposição QR de A usando os
    // Refletores de Householder.

    [m, n] = size(A)
    
    k = min(m, n)
    R = A
    U = zeros(m, n)

    // Construindo vetor de Househoulder
    for j = 1:k
        y = R(j:m,j)
        
        // Calculando matriz de Houlsehoulder para essa coluna.
        y_norm = norm(y)
        y_len = size(y, 'r')

        e_1 = zeros(y_len, 1)
        e_1(1) = 1

        if y(1) < 0 then
            y_norm = -y_norm
        end

        u = y - y_norm * e_1

        H_j = eye(y_len, y_len) - 2 * u * u'

        // Matriz H ampliada para aplicar em toda a matriz A.
        H = eye(m, m)
        H(k:m, k:m) = H_j

        // Atualizando R.
        R = H * R
    end
endfunction

function [Q] = Constroi_Q_House(U)
    // Constrói a matriz Q da decomposição QR a partr da matriz U do método
    // de Householder.

    [m, k] = size(Q)

    Q = eye(m, m)

    for i = 1:k
        u = U(:, i)
        Q = Q - 2 * Q * u * u'
    end
endfunction

A = [1 1 1;
     1 0 1;
     0 1 1]

[U, R] = qr_House(A)
[Q] = Constroi_Q_House(U)
disp(Q, R)

// Conferindo se funcionou.
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp(sum(Q' * Q))
disp("Q * R")
P_inv = inv(P)
disp(Q * R * P_inv) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)