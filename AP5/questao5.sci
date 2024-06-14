function [U, R] = qr_House_v2(A)
	// Description of qr_House(input)

	// Pegando dimensões de A.
	[m, n] = size(A)
	k = min(m - 1, n)

	// Inicializando matrizes U e R.
	U = zeros(m, k)
	R = A

	for j = 1:k
		// Selecionando coluna atual de A (R).
		x = R(j:m, j)
		x_size = size(x, 'r')
		
		// Criando vetor u.
		alpha = norm(x)
		if x(1) >=0 then
			alpha = -alpha
		end
		e_1 = [1; zeros(x_size - 1, 1)]
		u = x - alpha * e_1
		u = u/norm(u)
		// Adicionando vetor u à matriz U.
		U(j:m, j) = u

		// Aplicando reflexão.
		R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n))
		
		// R = H * R
	end

endfunction

function [Q] = Constroi_Q_House(U)
	// Cria matriz ortogonal Q a partir da matriz U.

	[m, n] = size(U)
	k = min(m, n)
	Q = eye(m, m)
	for j = n:-1:1
		u = U(j:m, j)
		u_size = size(u, 'r')

		// Aplicando reflexão.
		Q(j:m, :) = Q(j:m, :) - 2 * u * (u' * Q(j:m, :))
	end
endfunction

function [S] = espectro(A, tol)
    // Computes the eigenvalues of the matrix.
    
    // Número de linnhas de A.
    n = size(A, 'r')
    // Inicializando A_k.
    A_k = A
    // Inicializando vetores com autovalores.
    eig_1 = diag(A_k)
    eig_0 = eig_1

    tolerancia_alcancada = %f

    while ~tolerancia_alcancada then
        // Calculando Q e R.
        [U, R] = qr_House_v2(A_k)
        [Q] = Constroi_Q_House(U)
        
        A_k = R * Q
        eig_1 = diag(A_k)

        erro = norm((eig_1 - eig_0), 2)

        if erro < tol then
            tolerancia_alcancada = %t
        end

        eig_0 = eig_1
    end

    S = eig_1
endfunction

A = [1 1 0;
     1 1 1;
     0 1 1]

[S] = espectro(A, 0.00001)
disp("S")
disp(S)
disp("Espectro do scilab: ")
disp(spec(A))

// Matriz mal-condicionada.
A = [44 5 26 23;
     5 1 7 3;
     26 7 31 13;
     23 3 13 12]

[S] = espectro(A, 0.00001)
disp("S")
disp(S)
disp("Espectro do scilab: ")
disp(spec(A))