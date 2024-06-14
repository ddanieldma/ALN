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

        // Subtraindo a projeção para encontrar o perpendicular.
        v = v - projecao
    end

    v_norm = norm(v, 2)

    // Normalizando vetor.
    q = v/v_norm
    r_vetor = [r_vetor; v_norm]

endfunction

function [Q, R] = qr_GS(A)
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

// Matriz simples.
A = [1 1 1;
     1 0 1;
     0 1 1]

disp("A")
disp(A)
[Q, R] = qr_GS(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

// Matriz aleatória.
A = [16.90  5.25 13.27 22.56 10.08;
8.06 23.25  3.21 21.08 13.28;
20.65 10.43 19.83 11.52  5.40;
21.78  6.29 17.63 17.70 20.36;
23.74 10.50  1.78 10.18  4.68]
disp("A")
disp(A)

[Q, R] = qr_GS(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)