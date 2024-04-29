function [sol, final_it, k, norm_res] = GS_Method(A, b, x_0, E, M, norm_type)
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

A = [3, -2, 1;
     1, 3, 2;
     -1, 2, 4]

initial_vector = [1; 1; 1]
max_iterations = 20
tolerance = 0.01
b = [1; 1; 1]


// com a norma soma
// disp("Norma soma")
// norm_type = 1 // norma soma
// [x, final_it, num_it, norm_res] = GS_Method(A, b, initial_vector, tolerance, max_iterations, norm_type)

// // disp("x final")
// // disp(x)

// // disp("Norma da ultima iteração")
// // disp(final_it)

// disp("Número de iterações")
// disp(num_it)

// disp("Norma do resíduo")
// disp(norm_res)


// // com a norma euclidiana
// disp("Norma euclidiana")
// norm_type = 2 // norma euclidiana
// [x, final_it, num_it, norm_res] = GS_Method(A, b, initial_vector, tolerance, max_iterations, norm_type)

// disp("x final")
// disp(x)

// disp("Norma da ultima iteração")
// disp(final_it)

// disp("Número de iterações")
// disp(num_it)

// disp("Norma do resíduo")
// disp(norm_res)


// com a norma máximo
disp("Norma máximo")
norm_type = %inf // norma máximo
[x, final_it, num_it, norm_res] = GS_Method(A, b, initial_vector, tolerance, max_iterations, norm_type)

disp("x final")
disp(x)

disp("Norma da ultima iteração")
disp(final_it)

disp("Número de iterações")
disp(num_it)

disp("Norma do resíduo")
disp(norm_res)