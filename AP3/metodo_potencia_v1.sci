function [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, epsilon, M)
	// Calcula o autovalor dominante da matriz de maneira iterativa usando o método
	// da potência.

	k=0
	// Dividindo vetor pela valor da coordenada de maior módulo para que essa
	//seja igual a 1.
	x0 = x0/max(abs(x0))
	x1 = A*x0
	n_erro = epsilon + 1 // Obriga a entrar no loop

	while k<=M && n_erro >= epsilon then
		lambda = max(abs(x1))
		x1 = x1/lambda
		n_erro = norm((x1 - x0), %inf)
		x0 = x1
		x1 = A*x0
		k = k+1
	end
endfunction

// A = [1 1;
// 	 2 0]

// x0 = [1; 0]

// [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

// disp("Matriz A")
// disp(A)
// disp("Vetor inicial")
// disp(x0)
// disp("Autovalor dominante")
// disp(lambda)
// disp("Autovetor dominante")
// disp(x1)
// disp("Número de iterações")
// disp(k)
// disp("Erro das duas útlimas iterações")
// disp(n_erro)

// A = [1 2 3;
// 	 2 1 4;
// 	 3 4 1]

// x0 = [1; 0; 0]

// [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

// disp("Matriz A")
// disp(A)
// disp("Vetor inicial")
// disp(x0)
// disp("Autovalor dominante")
// disp(lambda)
// disp("Autovetor dominante")
// disp(x1)
// disp(x1/max(x1))
// disp("Número de iterações")
// disp(k)
// disp("Erro das duas útlimas iterações")
// disp(n_erro)

A = [-2.5 3.5 -2.5;
	 -1 2 1;
	 -3.5 3.5 -1.5]

x0 = [1; 0; 0]

[lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

disp("Matriz A")
disp(A)
disp("Vetor inicial")
disp(x0)
disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(x1)
disp(x1/max(x1))
disp("Número de iterações")
disp(k)
disp("Erro das duas útlimas iterações")
disp(n_erro)