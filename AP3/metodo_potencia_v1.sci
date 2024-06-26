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

// A = [1 1;
// 	 2 0]

// x0 = [1; 0]

// [lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

// disp("Autovalor dominante")
// disp(lambda)
// disp("Autovetor dominante")
// disp(Ajusta_vetor(x1))
// disp("Número de iterações")
// disp(k)
// disp("Erro da última iteração")
// disp(n_erro)
// disp("Autovalores pelo scilab")
// disp(spec(A))

// A = [1 2 3 4;
// 	 2 1 4 10;
// 	 3 4 1 -7;
// 	 5 4 5 1]

// x0 = [1; 0; 0; 0]

// for i=1:1:5
// 	disp("===================================")
	
// 	epsilon = 1/(10**i)
// 	disp("epsilon:")
// 	disp(epsilon)

// 	[lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, epsilon, 100)
	
// 	disp("Aproximação do autovalor")
// 	disp(lambda)
// 	disp("Autovetor dominante")
// 	disp(x1/max(x1))
// 	disp("Número de iterações")
// 	disp(k)
// end

// disp("Autvalor dominante da matriz")
// autovalores = spec(A)
// [dom, idx] = max(abs(autovalores))
// if autovalores(idx) < 0 then
// 	dom = -dom
// end
// disp(dom)

A = [-2.5 3.5 -2.5;
	 -1 2 1;
	 -3.5 3.5 -1.5]

x0 = [1; 0; 0]

disp("Autovalores da matriz A")
disp(spec(A))

[lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(Ajusta_vetor(x1))
disp("Número de iterações")
disp(k)
disp("Erro da última iteração")
disp(n_erro)

// Pra um autovalor dominante negativo o método parece não ter convergido,
// ultrapassando o limite de iterações mas encontrando o valor absoluto do autovalor
// e encontrando o autovetor dominante.
// Fazendo a correção de verificar se o valor dominante no vetor é negativo
// e transformar o lambda para o negativo caso seja funfou.