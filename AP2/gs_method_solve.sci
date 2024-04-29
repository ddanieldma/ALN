function [x] = L_solve(L, b)
	[m, n] = size(b)
    [p, q] = size(L)
    
    x=zeros(m, n);
    
    x(1)=b(1)/L(1,1)
    
    for i=2:m
        x(i)=(b(i)-L(i,1:(i-1))*x(1:(i-1)))/L(i,i);
    end
endfunction

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
		
		LD = L+D
		b_solve = -1*U*x_k+b
		[x_k1] = L_solve(LD, b_solve)
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

initial_vector = [0; 0; 0]
max_iterations = 20
tolerance = 0.01
b = [1; 1; 1]


// // com a norma soma
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

// // disp("x final")
// // disp(x)

// // disp("Norma da ultima iteração")
// // disp(final_it)

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

// mesma coisa de antes, a norma infinito parece que toma passos menores
// a cada iteração e com o mesmo número de iterações termina com um
// resíduo menor.

// parece ter pouca diferença no resultado da calculando com a matriz inversa
// ou resolvendo. Contudo, no geral parece que calcular com a matriz inversa
// demanda mais iterações.