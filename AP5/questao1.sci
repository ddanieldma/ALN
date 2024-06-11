function [q, r_vetor] = Project_into(Q, a)
    // Projeto um vetor no espaço definido pelos vetores dados.
    // q é o vetor no qual estamos projetando a.
    // Retorna a projeção e a coluna de r definida por ela.

    // Número de colunas.
    cols = size(Q, 'c')
    r_vetor = 0

    // Iniciando operação que vai definir v.
    v = a
    for i = 1:cols
        q = Q(:, i)
        coeficiente_projecao = q' * a
        projecao = coeficiente_projecao * q
        if i==1 then
            r_vetor = coeficiente_projecao
        else
            // Adicionando coeficiente na matriz de coeficientes.
            r_vetor = [r_vetor; coeficiente_projecao]
        end

        // Subtraindo a projeção para enconra
        v = v - projecao
    end

    v_norm = norm(v, 2)

    // Normalizando vetor.
    q = v/v_norm
    r_vetor = [r_vetor; v_norm]

endfunction

function [Q, R] = qr_GS(A)
    // Calcula a decomposição QR da matriz A.
    
    // Tamanho de A.
    n = size(A, 'r')

    // Inicializando Q e R.
    Q = zeros(n, n)
    R = zeros(n, n)

    for j = 1:n
        a_j = A(:, j)
        // Vetor v começa como a.
        v = a_j
        for i = 1:j-1
            // Pegando q's anteriores.
            q_i = Q(:, i)
            
            // Calculando coeficientes.
            r_ij = q_i' * a_j
            R(i, j) = r_ij

            // Calculando projeção de a no espaço.
            v = v - r_ij * q_i
        end
        v_norm = norm(v, 2)
        R(j, j) = v_norm
        
        // Normalizando v.
        Q(:,j) = v/v_norm
    end
endfunction

function [Q, R] = My_qr_GS(A)
    // Calcula a decomposição QR da matriz A.
    
    n = size(A, 'r');

    // Primeiro vetor é a primeira coluna de A.
    v_current = A(:, 1)
    // Normalizando.
    v_norm = norm(v_current, 2)
    q_current = v_current/v_norm
    r_current = [v_norm; zeros(n - 1, 1)]

    Q = q_current
    R = [r_current]

    for i = 2:n
        // Os outros vetores são a projeção da próxima coluna de A no espaço
        // formado.

        // Calculando proxima projeção.
        [q_current, r_vetor] = Project_into(Q, A(:, i))

        // Atualizando variáveis
        Q = [Q q_current]
        // Adicionando zeros em r_vetor caso seja necessário.
        r_vetor = [r_vetor; zeros(n-i,1)]
        R = [R r_vetor]
    end
endfunction

function [error] = Compute_accuracy_QR(Q, R, A)
    // Compute how close to the initial A matrix QR is.

    QR = Q * R
    error = A - QR
    error = sum(error)
endfunction

// A = [1 1 1;
//      1 0 1;
//      0 1 1]

// [Q, R] = My_qr_GS(A)

// disp(Q, R)

// // Conferindo se funcionou.
// disp(Q * R) // Tem que printar A.
// disp(Q' * Q) // Tem que printar a identidade.
// disp("Erro: ")
// erro = Compute_accuracy_QR(Q, R, A)
// disp(erro)

// A = [1 1 1;
//      1 0 1;
//      0 1 1]

// [Q, R] = qr_GS(A)

// disp(Q, R)

// // Conferindo se funcionou.
// disp(Q * R) // Tem que printar A.
// disp(Q' * Q) // Tem que printar a identidade.
// disp("Erro: ")
// erro = Compute_accuracy_QR(Q, R, A)
// disp(erro)
// Erro suficientemente pequeno para ser ignorado.

A = [1 10 8 2;
     5 6 3 7;
     4 2 2 0;
     2 0 8 5]

[Q, R] = qr_GS(A)

disp(Q, R)

// Conferindo se funcionou.
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

A = [16.90  5.25 13.27 22.56 10.08;
 8.06 23.25  3.21 21.08 13.28;
20.65 10.43 19.83 11.52  5.40;
21.78  6.29 17.63 17.70 20.36;
23.74 10.50  1.78 10.18  4.68]

[Q, R] = qr_GS(A)

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