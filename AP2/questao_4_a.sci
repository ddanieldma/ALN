function [sol, final_it, k, norm_res] = Jacobi_Method(A, b, x_0, E, M, norm_type)

    [n]=size(A,1)

    L = tril(A, -1)
    U = triu(A, 1)
    
    D = zeros(n, n)
    
    elementos_diagonal = diag(A)

    for i=1:n
        D(i,i) = elementos_diagonal(i)
    end

    k=0
    D_inv = inv(D)
    M_j = (-D_inv)*(L+U)
	disp("raio espectral")
	disp(max(abs(spec(M_j))))
    c_j = D_inv*b
    x_k = x_0

    while k < M
        k = k + 1
        x_k1 = M_j*x_k + c_j
        final_it = norm((x_k1 - x_k), norm_type)
        x_k = x_k1

        if final_it < E then
            disp("Passou a tolerância")
            break
        end
        
        if k == M then
            disp("Ultrapassou o número máximo de iterações")
        end
    end

    sol = x_k1

    res = norm((b - A*x_k1), norm_type)
    norm_res = norm(res, norm_type)

endfunction

A = [2 -1 1;
	 2 2 2;
	 -1 -1 2
]
b = [-1; 4; -5]

x_0 = [0; 0; 0]
E = 0.001
M = 26
norm_type = 2

[aprox, final_it, k, norm_res] = Jacobi_Method(A, b, x_0, E, M)

disp("Aproximação")
disp(aprox)

disp("Norma do resíduo")
disp(norm_res)

// a matriz não tem diagonal dominante e nem raio espectral menor que 1.
// Como o raio não é muito maior que 1, os números não estouram tanto.