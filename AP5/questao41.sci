
function [U, R] = qr_House_v1(A)
    // Decompõe matriz A em U e R usando refletores de Househoulder.
    
    // Pegando dimensões de A.
    [m, n] = size(A)
    R = A
    U = zeros(m, n)

    for j = 1:n
        x = R(j:m, j)
        x_size = size(x, 'r')
		
		// Criando vetor u.
		alpha = norm(x)
		if x(1) >=0 then
			alpha = -alpha
		end
		e_1 = [1; zeros(x_size - 1, 1)]
		u = x - alpha * e_1
        u_norm = norm(u)

        // Evitando Nan.
        if u_norm > 0 then
            // Normalizando U.
            u = u/u_norm
            // Aplicando reflexão.
            R(j:m, j:n) = R(j:m, j:n) - 2 * u * (u' * R(j:m, j:n))
            U(j:m, j) = u
        end
    end
endfunction

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

// Matriz simples.
A = [1 1 1;
     1 0 1;
     0 1 1]
// 
disp("Primeira versão")
[U, R] = qr_House_v1(A)
[Q] = Constroi_Q_House(U)

disp(U, Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

disp("Segunda versão")
[U, R] = qr_House_v1(A)
[Q] = Constroi_Q_House(U)

disp(U, Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

// Matriz aleatória.
A = [16.90  5.25 13.27 22.56 10.08;
8.06 23.25  3.21 21.08 13.28;
20.65 10.43 19.83 11.52  5.40;
21.78  6.29 17.63 17.70 20.36;
23.74 10.50  1.78 10.18  4.68]

disp("Primeira versão")
[U, R] = qr_House_v1(A)
[Q] = Constroi_Q_House(U)

disp(U, Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

disp("Segunda versão")
[U, R] = qr_House_v1(A)
[Q] = Constroi_Q_House(U)

disp(U, Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)