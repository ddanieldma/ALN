// Função do método iterativo de jacobi para aproximar x em Ax = b.
function [sol, final_it, k, norm_res] = Jacobi_Method(A, b, x_0, E, M, norm_type)
    // Entradas:
    // A: matriz de coeficientes
    // b: vetor de respostas
    // x_0: aproximação incial
    // E: tolerância
    // M: número máximo de iterações
    // norm_type: tipo de norma a ser utilizado (1, 2, ou %inf)

    // Saídas:
    // sol: A solução x_k encontrada
    // final_it: A norma da diferença entre as duas últimas aproximações
    // k: O número de iterações efetuadas
    // res: A norma do resíduo

    // Critério de parada: norma da diferença entre as duas útlimas
    // aproximações menor que E ou número de iterações maior que M

    [n]=size(A,1)

    // Criando matrizes L,D e U.
    L = tril(A, -1)
    U = triu(A, 1)
    
    D = zeros(n, n)
    
    elementos_diagonal = diag(A)

    for i=1:n
        D(i,i) = elementos_diagonal(i)
    end

    sol = 0
    res = 0

    D_inv = inv(D)

    x_k = x_0
    while k < M
        k = k + 1
        x_k1 = -D * (L + U)x_k + D_inv * B
        final_it = norm((x_k1 - x_k), norm_type)

        if final_it < E then
            break
        end
    end

    sol = x_k1

    res = norm((b - A*x_k1), norm_type)
    norm_res = norm(res, norm_type)

endfunction

A = [1, 2, 3;
     4, 5, 6;
     7, 8, 9]

[a, b] = Jacobi_Method(A, 1, 2, 3, 4)
