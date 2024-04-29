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
	// matriz do método
	M_gs = (-LD_inv)*U
	// constante do método
	c_gs = LD_inv*b
    x_k = x_0

    while k < M
        k = k + 1
		
		x_k1 = M_gs*x_k + c_gs
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

function [A, b, x_0] = gerar_matrizes(n)
    A = rand(n, n);
    x = rand(n, 1);
    x_0 = zeros(n, 1);

    for i = 1:n
        for j = 1:n
            A(i, i) = A(i,i) + A(i, j) + A(j, i)
        end
    end

    b = A*x
endfunction

tolerance = 0.001
max_it = 1000
norm_type = 2

n = 10
[A, b, x_0] = gerar_matrizes(n)

disp("========================================================")
disp("n = ", n)
tic()
GS_Method_Std(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS normal")
disp(t)

tic()
GS_Method_Solve(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS com eliminação")
disp(t)


n = 100
[A, b, x_0] = gerar_matrizes(n)

disp("========================================================")
disp("n = ", n)
tic()
GS_Method_Std(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS normal")
disp(t)

tic()
GS_Method_Solve(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS com eliminação")
disp(t)

n = 1000
[A, b, x_0] = gerar_matrizes(n)

disp("========================================================")
disp("n = ", n)
tic()
GS_Method_Std(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS normal")
disp(t)

tic()
GS_Method_Solve(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS com eliminação")
disp(t)

n = 3000
[A, b, x_0] = gerar_matrizes(n)

disp("========================================================")
disp("n = ", n)
tic()
GS_Method_Std(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS normal")
disp(t)


tic()
GS_Method_Solve(A, b, x_0, tolerance, max_it, norm_type)
t = toc()
disp("Tempo GS com eliminação")
disp(t)
