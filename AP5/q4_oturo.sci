function [U, R] = qr_House_v1(A)
    // Inicializa as dimensões da matriz A
    [m, n] = size(A);
    R = A;
    U = zeros(m, n);

    // Loop sobre cada coluna de A
    for k = 1:n
        x = R(k:m, k);
        e = zeros(length(x), 1);
        e(1) = norm(x);
        
        // Calcula o vetor u
        u = x - e;
        norm_u = norm(u);
        // Verifica se norm_u é diferente de zero para evitar NaN
        if norm_u > 0 then
            u = u / norm_u;
            // Aplica a reflexão de Householder
            R(k:m, k:n) = R(k:m, k:n) - 2 * u * (u' * R(k:m, k:n));
            U(k:m, k) = u;
        end
    end
end

function [U, R] = qr_House_v2(A)
    // Inicializa as dimensões da matriz A
    [m, n] = size(A);
    k = min(m-1, n);
    R = A;
    U = zeros(m, k);

    // Loop sobre cada coluna de A
    for j = 1:k
        x = R(j:m, j);
        e = zeros(length(x), 1);
        e(1) = norm(x);
        u = x - e;
        u = u / norm(u);

        // Aplica a reflexão de Householder
        R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n));
        U(j:m, j) = u;
    end
end

function [Q] = constroi_Q_House_v1(U)
    // Inicializa as dimensões da matriz U
    [m, n] = size(U);
    Q = eye(m,m);

    // Loop inverso para construir a matriz Q
    for k = n:-1:1
        u = U(k:m, k);
        // Verifica se u não é um vetor nulo
        if norm(u) > 0 then
            Q(k:m, :) = Q(k:m, :) - 2 * u * (u' * Q(k:m, :));
        end
    end
end

function [Q] = constroi_Q_House_v2(U)
    // Inicializa as dimensões da matriz U
    [m, k] = size(U);
    Q = eye(m,m);

    // Loop inverso para construir a matriz Q
    for j = k:-1:1
        u = U(j:m, j);
        Q(j:m, :) = Q(j:m, :) - 2 * u * (u' * Q(j:m, :));
    end
end

function [U, R] = my_qr_House_v2(A)
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

		// Criando matriz de Househoulder da iteração atual.
		// H_j = eye(x_size, x_size) - 2 * u * u'
		// Aplicando reflexão.
		R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n))
		
		// R = H * R
	end

endfunction

function [Q] = my_Constroi_Q_House_v2(U)
	// Description of Constroi_QR(U)
	[m, n] = size(U)
	k = min(m, n)
	Q = eye(m, m)
	for j = n:-1:1
		u = U(j:m, j)
		u_size = size(u, 'r')

		// Aplicando reflexão.
		Q(j:m, :) = Q(j:m, :) - 2 * u * (u' * Q(j:m, :))

		// H_j = eye(u_size, u_size) - 2 * u * u'
		// H = eye(k, k)
		// H(j:k, j:k) = H_j
		// Q = H * Q
	end
endfunction


A = [1 1 1 0;
     1 0 1 1;
     0 1 1 1]

[U, R] = my_qr_House_v2(A)
disp("U")
disp(U)
[Q] = my_Constroi_Q_House_v2(U)

disp("Q * R")
disp(Q * R)

disp("QT * Q")
disp(Q' * Q)