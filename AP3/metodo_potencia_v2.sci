function output = Metodo_potencia_v2(input)
	// Encontra o autovalor dominante da matriz pelo método da potência, utilizando
	// o Quociente de Rayleigh.

	k = 0
	x0 = x0/norm(x0, 2)
	x1 = A*x0
	n_erro = epsilon+1

	while k<=M && n_erro >= epsilon then
		lambda = x1*x0
		if lambda<0 then
			x1 = -x1
		end
		x1 = x1/norm(x1, 2)
		n_erro = norm((x1-x0), 2)
		x0 = x1
		x1 = A*x0
		k = k+1
	end

endfunction

A = [1 1;
	 2 0]

x0 = [1; 0]

[lambda, x1, k, n_erro] = Metodo_potencia_v1(A, x0, 0.001, 100)

disp("Matriz A")
disp(A)
disp("Vetor inicial")
disp(x0)
disp("Autovalor dominante")
disp(lambda)
disp("Autovetor dominante")
disp(x1)
disp("Número de iterações")
disp(k)
disp("Erro das duas útlimas iterações")
disp(n_erro)