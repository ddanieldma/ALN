// function [projecao] = Project_into(q, a)
//     // Projeto um vetor no espaço definido pelos vetores dados.
//     // q é o vetor no qual estamos projetando a.

//     numerador = q * a
//     denominador = q * q

//     projecao = numerador/denominador * q

// endfunction

function [q, r_vetor] = Other_project_into(Q, a)
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

// function [Q, R] = ComputeQR(A)
//     // Calcula a decomposição QR da matriz A.
    
//     n = size(A, 'c');

//     // Primeiro vetor é a primeira coluna de A.
//     v_current = A(:,1)
//     q_current = v_current/norm(v_current, 2)
//     v_previous = v_current
//     q_previous = q_current

//     Q = q_1
//     R = []

//     for i = 2:n
//         // Os outros vetores são a projeção da próxima coluna de A no espaço
//         // formado.

//         // Calculando proxima projeção.
//         v_current = v_previous - proj(q_previous)
//         q_current = v_current/norm(v_current, 2)
//         Q = [Q, q_current]

//         // Atualizando variáveis
//         v_previous = v_current
//         q_previous = q_current
//     end
// endfunction

function [Q, R] = Other_compute_QR(A)
    // Calcula a decomposição QR da matriz A.
    
    n = size(A, 'c');

    // Primeiro vetor é a primeira coluna de A.
    v_current = A(:,1)
    // Normalizando.
    v_norm = norm(v_current, 2)
    q_current = v_current/v_norm
    r_current = [v_norm; zeros(n-1, 1)]
    
    // Atualizando variáveis.
    v_previous = v_current
    q_previous = q_current

    Q = q_current
    R = [r_current]

    disp(R)

    for i = 2:n
        // Os outros vetores são a projeção da próxima coluna de A no espaço
        // formado.

        // Calculando proxima projeção.
        [q_current, r_vetor] = Other_project_into(Q, A(:, i))
        if i==2 then
            disp(r_vetor)
        end

        // Atualizando variáveis
        Q = [Q q_current]
        // Adicionando zeros em r_vetor caso seja necessário.
        r_vetor = [r_vetor; zeros(n-i,1)]
        R = [R r_vetor]
    end
endfunction

A = [1 1 1;
     1 0 1;
     0 1 1]

[Q, R] = Other_compute_QR(A)

disp(Q, R)

// Conferindo se funcionou.
disp(Q * R) // Tem que printar A.
disp(Q' * Q) // Tem que printar a identidade.