function [Q, R] = qr_GSM(A)
	// Outra forma de calcular a decomposição QR da matriz A.

	n = size(A, 'r')

	for j = 1:n
		v = A(:, j)
		for i = 1:j-1
			q_i = Q(:, i)
			r_ij = q_i' * v
			R(i, j) = r_ij
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
    error = sum(error)
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

// Matriz aleatória.
A = [16.90  5.25 13.27 22.56 10.08;
8.06 23.25  3.21 21.08 13.28;
20.65 10.43 19.83 11.52  5.40;
21.78  6.29 17.63 17.70 20.36;
23.74 10.50  1.78 10.18  4.68]

[Q, R] = qr_GSM(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A)
disp(erro)