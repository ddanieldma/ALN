
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
	// Decompõe matriz A em U e R usando refletores de Househoulder.

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

function [error] = Compute_accuracy_QR(Q, R, A)
    // Compute how close to the initial A matrix QR is.

    QR = Q * R
    error = A - QR
    error = sum(abs(error))
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
[U, R] = qr_House_v2(A)
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
// Matriz aleatória com colunas quase lineramente dependentes.
A = [ 51. 92. 14. 71. 143.0000242 ;
     20. 82. 86. 74. 101.99980867;
     87. 99. 23.  2. 185.99982751;
     52.  1. 87. 29. 52.99994377;
      1. 63. 59. 20. 63.99989872]
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
[U, R] = qr_House_v2(A)
[Q] = Constroi_Q_House(U)

disp(U, Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)