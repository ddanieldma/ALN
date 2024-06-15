function [P] = Create_perm_col(size_m, col_1, col_2)
    //create identity matrix
    I = eye(size_m, size_m)
    // permutate its lines
    I(:, [col_1, col_2]) = I(:, [col_2, col_1])
    P = I
endfunction

// Decomposição QR com pivoteamento de colunas.
function [Q, R, P] = qr_GSP(A)
    // Calcula a decomposição QR da matriz A usando pivoteamento de colunas.
    
    // Tamanho de A.
    n = size(A, 'r')

    // Inicializando Q e R.
    Q = zeros(n, n)
    R = zeros(n, n)

	// Inicializando matriz que vai armazenar as permutações.
	P = eye(n, n)

    for j = 1:n
        // Fazendo pivoteamento.
		// Calculando normas.
		normas = zeros(1, n)
		for i = 1:n
			normas(i) = norm(A(:,i))
		end
		// Obtendo maior coluna e seu índice.
		[_, coluna_indice] = max(normas);
        coluna_indice = coluna_indice + j - 1
		// Criando permutação.
        [Perm] = Create_perm_col(n, j, coluna_indice)
		// Permutando A.
		A = A * Perm
		// Armezando permutação.
		P = P * Perm

		a_j = A(:, j)
        // Vetor v começa como a.
        v = a_j
        for i = 1:j-1			
			// Pegando q's anteriores.
            q_i = Q(:, i)
            
            // Calculando coeficientes.
            r_ij = q_i' * a_j
            R(i, j) = r_ij

            // Calculando projeção de a no espaço.
            v = v - r_ij * q_i
        end
        v_norm = norm(v)
        R(j, j) = v_norm
        
        // Normalizando v.
        Q(:,j) = v/v_norm
    end
endfunction

function [error] = Compute_accuracy_QR(Q, R, A, P)
    // Compute how close to the initial A matrix QR is.

    QR = Q * R
    QR = QR * P
    error = A - QR
    error = sum(abs(error))
endfunction

// Matriz simples.
A = [1 1 1;
     1 0 1;
     0 1 1]

[Q, R, P] = qr_GSP(A)

disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R * P) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A, P)
disp(erro)

// Matriz aleatória com colunas quase lineramente dependentes.
A = [ 51. 92. 14. 71. 143.0000242 ;
     20. 82. 86. 74. 101.99980867;
     87. 99. 23.  2. 185.99982751;
     52.  1. 87. 29. 52.99994377;
      1. 63. 59. 20. 63.99989872]
//   

[Q, R, P] = qr_GSP(A)
disp(Q, R)

disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R * P) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, A, P)
disp(erro)