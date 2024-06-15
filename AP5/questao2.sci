function [Q, R] = qr_GSM(A)
	// Outra forma de calcular a decomposição QR da matriz A.

	// Pegando número de linhas de A.
	n = size(A, 'r')

	for j = 1:n
		// Vetor atual é a j-esima coluna de A.
		v = A(:, j)
		
		// A partir da segunda iteração:
		for i = 1:j-1
			// Fazendo projeções de v no espaço dos q's já criados.
			q_i = Q(:, i)
			// Coeficiente.
			r_ij = q_i' * v
			R(i, j) = r_ij
			// Pegando perpendicular da projeção.
			v = v - r_ij * q_i
		end

		v_norm = norm(v, 2)
		// Normalizando v;
		Q(:, j) = v/v_norm
		R(j, j) = v_norm
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

[Q, R] = qr_GSM(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)

// Matriz aleatória com colunas quase lineramente dependentes.
A = [ 51. 92. 14. 71. 143.0000242 ;
     20. 82. 86. 74. 101.99980867;
     87. 99. 23.  2. 185.99982751;
     52.  1. 87. 29. 52.99994377;
      1. 63. 59. 20. 63.99989872]
//   

[Q, R] = qr_GSM(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)