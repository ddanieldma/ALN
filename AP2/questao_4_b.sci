// b)
function [sol, final_it, k, norm_res] = GS_Method(A, b, x_0, E, M, norm_type)
    [n]=size(A,1)

    L = tril(A, -1)
    U = triu(A, 1)
    
    D = zeros(n, n)
    
    elementos_diagonal = diag(A)

    for i=1:n
        D(i,i) = elementos_diagonal(i)
    end

    k=0

	LD_inv = inv(L + D)
	M_gs = (-LD_inv)*U
	c_gs = LD_inv*b
    x_k = x_0

    while k < M
        k = k + 1
		
		x_k1 = M_gs*x + c_gs
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
// matriz não tem diagonal estritamente dominante

// conferindo se A é positiva definida
disp("Autovalores")
disp(spec(A))
// parece ser

b = [-1; 4; -5]
x_0 = [0; 0; 0]
E = 0.00001
M = 100000
norm_type = %inf

[aprox, final_it, k, norm_res] = GS_Method(A, b, x_0, E, M)

disp("Aproximação")
disp(aprox)

disp("Norma do resíduo")
disp(norm_res)

disp("Iterações")
disp(k)

// é esquisito que foram poucas iterações e a norma resíduo parece muito 
// grande