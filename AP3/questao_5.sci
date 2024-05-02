function [x] = Ajusta_vetor(x)
	// Retorna o vetor ajustado, dividindo-o pelo seu elemento da maior
	// valor absoluto.

	[max_value, max_index] = max(abs(x))
	if x(max_index) < 0 then
		max_value = -max_value
	end
	x = x/max_value
endfunction

function [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, epsilon, M)
	// Calcula o autovalor dominante da matriz de maneira iterativa usando o método
	// da potência.

	// Inicializando contador de iterações.
	k=0
	
	// Dividindo vetor pela valor da coordenada de maior módulo para que essa
	//seja igual a 1.
	[max_value, max_index] = max(abs(x0))
	if x0(max_index) < 0 then //
		lambda = -lambda //
	end //

	x0 = x0/max_value
	// Codigo adicionado pra corrigir erro.
	
	// Primeira iteração.
	x1 = A*x0
	
	// Iniciando erro maior que epsilon para que o codigo entre no loop.
	n_erro = epsilon + 1 // Obriga a entrar no loop

	while k<=M && n_erro >= epsilon then
		// A aproximação do autovalor é o elemento de maior valor absoluto
		// do vetor.

		// Apesar de precisarmos do elemento com maior valor absoluto
		// a operação não deve ser feita com o valor absoluto.
		[lambda, max_index] = max(abs(x1)) //
		if x1(max_index) < 0 then //
			lambda = -lambda //
		end //

		// Redução de escala.
		x1 = x1/lambda

		// Erro da iteração atual.
		n_erro = norm((x1 - x0), %inf)
		
		// Indo para a próximo iteração
		x0 = x1
		x1 = A*x0
		k = k+1
	end
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


// Matriz com autovalores que são um o oposto do outro.

// A = [-2 -3 -3;
// 	 -4 -2 -4;
// 	 4 3 5]
// // 

// disp(spec(A))

// disp(spec(A))

// [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

// disp("Autovalor dominante")
// disp(lambda)
// disp("Autovetor dominante")
// disp(Ajusta_vetor(x1))
// disp("Número de iterações")
// disp(k)
// disp("Erro da última iteração")
// disp(n_erro)
// O método não converge, devido à inexistência de um autovalor dominante.

// Matriz com dois autovalores muito próximos.
A = [499/100 399/100 399/100;
	 -1/100 499/100 -1/100;
	 1/100 -399/100 101/100]
// 

x0 = [1; 0; 0]

[lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(Ajusta_vetor(x1))
disp("Número de iterações")
disp(k)

[lambda, x1, k, n_erro] = Metodo_potencia_v2(A, x0, 0.001, 100)

disp("=========================================")
disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(Ajusta_vetor(x1))
disp("Número de iterações")
disp(k)

// O algoritmo falha em convergir para algum dos autovalores e o autovetor final é uma combinação 
// dos dois autovetores.
// A versão 1 do método da potência não converge, apesar de acabar retornando um dos autovalores.
// A versão 2 converge para um autovalor "no meio dos dois" e retorna um autovetor que parece
// ser uma combinação dos dois autovetores, apesar de tender mais para um autovetor do que para o outro