// function [U, R] = qr_House(A)
//     // Matriz que a calcula a decomposição QR de A usando os
//     // Refletores de Householder.

//     [m, n] = size(A)
    
//     k = min(m, n)
//     R = zeros(m, n)
//     U = zeros(m, n)

//     for j = 1:k
//         x = A(j:m, j)
//         x_norm = norm(x, 2)
//         if x(1) < 0 then
//             x(1) = x(1) - x_norm
//         else
//             x(1) = x(1) + x_norm
//         end

//         u = x/x_norm
//         U(j:m, j) = u
//         A(j:m, j:n) = A(j:m, j:n) - 2 * u * (u' * A(j:m, j:n))
//     end

//     R = triu(A)
// endfunction

function [U, R] = qr_House(A)
    // Matriz que a calcula a decomposição QR de A usando os
    // Refletores de Householder.

    [m, n] = size(A)
    
    k = min(m, n)
    R = A
    U = zeros(m, n)

    // Construindo vetor de Househoulder
    for j = 1:k
        x = R(j:m, j)
        
        // Calculando matriz de Houlsehoulder para essa coluna.
        x_norm = norm(x)
        x_len = size(x, 'r')

        e_1 = zeros(x_len, 1)
        e_1(1) = 1

        if x(1) >= 0 then
            x_norm = -x_norm
        end

        u = x - x_norm * e_1
        u = u/norm(u)
        // Armazena u.
        U(j:m, j) = u

        H_j = eye(x_len, x_len) - 2 * u * u'

        // Matriz H ampliada para aplicar em toda a matriz R.
        H = eye(m, m)
        H(j:m, j:m) = H_j

        // Atualizando R.
        R = H * R
    end
endfunction

function [Q] = Constroi_Q_House(U)
    // Constrói a matriz Q da decomposição QR a partr da matriz U do método
    // de Householder.

    [m, n] = size(U)

    Q = eye(m, m)

    for i = n:-1:1
        u = U(i:m, i)
        H = eye(m, m)
        H(i:m, i:m) = eye(m - i + 1) - 2 * u * u'
        Q = Q * H
    end
endfunction

function [error] = Compute_accuracy_QR(Q, R, A)
    // Compute how close to the initial A matrix QR is.

    QR = Q * R
    error = A - QR
    error = sum(error)
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
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)