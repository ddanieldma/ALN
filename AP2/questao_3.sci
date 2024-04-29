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

function [sol, final_it, k, norm_res] = GS_Method_Std(A, b, x_0, E, M, norm_type)
    [n]=size(A,1)

    // Criando matrizes L,D e U.
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

function [x] = L_solve(L, b)
	[m, n] = size(b)
    [p, q] = size(L)
    
    x=zeros(m, n);
    
    x(1)=b(1)/L(1,1)
    
    for i=2:m
        x(i)=(b(i)-L(i,1:(i-1))*x(1:(i-1)))/L(i,i);
    end
endfunction

function [sol, final_it, k, norm_res] = GS_Method_Solve(A, b, x_0, E, M, norm_type)

    [n]=size(A,1)

    // Criando matrizes L,D e U.
    L = tril(A, -1)
    U = triu(A, 1)
    
    D = zeros(n, n)
    
    elementos_diagonal = diag(A)

    for i=1:n
        D(i,i) = elementos_diagonal(i)
    end

    k=0

	LD_inv = inv(L + D)
	// matriz do método
	M_gs = (-LD_inv)*U
	// constante do método
	c_gs = LD_inv*b
    x_k = x_0

    while k < M
        k = k + 1
		
		LD = L+D
		b_solve = -1*U*x_k+b
		[x_k1] = L_solve(LD, b_solve)
		// x_k1 = M_gs*x + c_gs
        final_it = norm((x_k1 - x_k), norm_type)
        x_k = x_k1

        // Parando o algoritmo quando a tolerância é passada
        if final_it < E then
            disp("Passou a tolerância")
            break
        end
        
        // Parando o algorítmo quando o número máximo de iterações
        // é ultrapassado
        if k == M then
            disp("Ultrapassou o número máximo de iterações")
        end
    end

    sol = x_k1

    res = norm((b - A*x_k1), norm_type)
    norm_res = norm(res, norm_type)

endfunction

A = [1 -4 2;
	 0 2 4;
	 6 -1 -2
]
// não tem raio espectral menor que 1 mas é estritamente dominante na
// diagonal.

b = [2; 1; 1]
x_0 = [0; 0; 0]
E = 0.001
M = 1000
norm_type = %inf

// não vai funcionar porque o raio espectral é maior que 1 e sem trocar
// linhas a matriz não tem diagonal estritametne dominante
disp("Jacobi: ")
[sol, final_it, k, norm_res] = Jacobi_Method(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)

disp("GS normal: ")
[sol, final_it, k, norm_res] = GS_Method_Std(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)

// provavelmente não resolve porque a eliminação com L não funciona sem troca
// de linhas

disp("GS resolvendo para L + D: ")
[sol, final_it, k, norm_res] = GS_Method_Solve(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)

// Antes de mudar a ordem das equações: apenas o método de Gauss Seidel funcionou.

// Mudando a ordem das equações para fazer de A diagonal
// dominante
aug = [A, b]

// trocando primeira e terceira linhas
aug([1, 3], :) = aug([3, 1], :)
// trocando segunda e terceira linhas
aug([2, 3], :) = aug([3, 2], :)

A = aug(:,1:3)
b = aug(:,4)

disp("Jacobi: ")
[sol, final_it, k, norm_res] = Jacobi_Method(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)
disp("GS normal: ")
[sol, final_it, k, norm_res] = GS_Method_Std(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)

disp("GS resolvendo para L + D: ")
[sol, final_it, k, norm_res] = GS_Method_Solve(A, b, x_0, E, M, norm_type) 

disp("X final")
disp(sol)

disp("Norma do resíduo")
disp(norm_res)

// antes de permutar as linhas apenas GS com a função de calcular inversa
// funcionou. Depois, todas funcionaram.