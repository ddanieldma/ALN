function [x] = Ajusta_vetor(x)
	// Retorna o vetor ajustado, dividindo-o pelo seu elemento da maior
	// valor absoluto.

	[max_value, max_index] = max(abs(x))
	if x(max_index) < 0 then
		max_value = -max_value
	end

	x = x/max_value
endfunction

function [lambda, x1, k, n_erro] = Metodo_potencia_v2(A, x0, epsilon, M)
	// Encontra o autovalor dominante da matriz pelo método da potência, utilizando
	// o Quociente de Rayleigh.

	// Inicializando contador de iterações.
	k = 0

	// Redução de escala.
	x0 = x0/norm(x0, 2)

	// Primeira iteração.
	x1 = A*x0

	// Inicializando erro.
	n_erro = epsilon+1

	while k<=M && n_erro >= epsilon then
		// Aproximando lambda com o quociente de Rayleigh.
		lambda = x1'*x0
		
		// Ajuste do vetor.
		if lambda<0 then
			x1 = -x1
		end

		// Redução de escala.
		x1 = x1/norm(x1, 2)

		// Erro da iteração atual.
		n_erro = norm((x1-x0), 2)

		// Indo para a próxima iteração.
		x0 = x1
		x1 = A*x0
		k = k+1
	end

endfunction

// A = [1 1;
// 	 2 0]

// x0 = [1; 0]

// [lambda, x1, k, n_erro] = Metodo_potencia_v2(A, x0, 0.001, 100)

// disp("Autovalor dominante")
// disp(lambda)
// disp("Autovetor dominante")
// disp(Ajusta_vetor(x1))
// disp("Número de iterações")
// disp(k)
// disp("Erro da última iteração")
// disp(n_erro)

// for i=1:1:5
// 	disp("===================================")
	
// 	epsilon = 1/(10**i)
// 	disp("epsilon:")
// 	disp(epsilon)

// 	[lambda, x1, k, n_erro] = Metodo_potencia_v2(A, x0, epsilon, 100)
	
// 	disp("Aproximação do autovalor")
// 	disp(lambda)
// 	disp("Autovetor dominante")
// 	disp(x1/max(x1))
// 	disp("Número de iterações")
// 	disp(k)
// end

A = [-2.5 3.5 -2.5;
	 -1 2 1;
	 -3.5 3.5 -1.5]

x0 = [1; 0; 0]

[lambda, x1, k, n_erro] = Metodo_potencia_v2(A, x0, 0.001, 100)

disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(Ajusta_vetor(x1))
disp("Número de iterações")
disp(k)
disp("Erro da última iteração")
disp(n_erro)