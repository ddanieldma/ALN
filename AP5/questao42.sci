// Algoritmos questao 1.
function [q, r_vetor] = Project_into(Q, a)
    // Projeto um vetor no espaço definido pelos vetores dados.
    // q é o vetor no qual estamos projetando a.
    // Retorna a projeção e a coluna de r definida por ela.

    // Número de colunas.
    cols = size(Q, 'c')
    r_vetor = 0

    // Iniciando operação que vai definir v.
    v = a
    for i = 1:cols
        q = Q(:, i)
        coeficiente_projecao = q' * a
        projecao = coeficiente_projecao * q
        if i==1 then
            r_vetor = coeficiente_projecao
        else
            // Adicionando coeficiente na matriz de coeficientes.
            r_vetor = [r_vetor; coeficiente_projecao]
        end

        // Subtraindo a projeção para encontrar o perpendicular.
        v = v - projecao
    end

    v_norm = norm(v, 2)

    // Normalizando vetor.
    q = v/v_norm
    r_vetor = [r_vetor; v_norm]

endfunction

function [Q, R] = qr_GS(A)
    // Calcula a decomposição QR da matriz A.
    
    n = size(A, 'r');

    // Primeiro vetor é a primeira coluna de A.
    v_current = A(:, 1)
    // Normalizando.
    v_norm = norm(v_current, 2)
    q_current = v_current/v_norm
    r_current = [v_norm; zeros(n - 1, 1)]

    Q = q_current
    R = [r_current]

    for i = 2:n
        // Os outros vetores são a projeção da próxima coluna de A no espaço
        // formado.

        // Calculando proxima projeção.
        [q_current, r_vetor] = Project_into(Q, A(:, i))

        // Atualizando variáveis
        Q = [Q q_current]
        // Adicionando zeros em r_vetor caso seja necessário.
        r_vetor = [r_vetor; zeros(n-i,1)]
        R = [R r_vetor]
    end
endfunction

function [error] = Compute_accuracy_QR(Q, R, A)
    // Compute how close to the initial A matrix QR is.

    QR = Q * R
    error = A - QR
    error = sum(error)
endfunction

// Algoritmos questao 2.
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

// Algoritmos questao 3.
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

// Algoritmos questao 4.

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

M1 = testmatrix('magi', 7)
H = testmatrix('hilb', 7)
M2 = testmatrix('magi', 6)

// --------------------------------------------------------
disp("Algoritmos questao 1")
disp("M1")
[Q, R] = qr_GS(M1)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M1)
disp(erro)

disp("M2")
[Q, R] = qr_GS(M2)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M2)
disp(erro)

disp("H")
[Q, R] = qr_GS(H)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, H)
disp(erro)

// --------------------------------------------------------
disp("Algoritmos questao 2")
disp("M1")
[Q, R] = qr_GSM(M1)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M1)
disp(erro)

disp("M2")
[Q, R] = qr_GSM(M2)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M2)
disp(erro)

disp("H")
[Q, R] = qr_GSM(H)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, H)
disp(erro)

// --------------------------------------------------------
disp("Algoritmos questao 3")
disp("M1")
[Q, R] = qr_GSP(M1)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M1)
disp(erro)

disp("M2")
[Q, R] = qr_GSP(M2)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M2)
disp(erro)

disp("H")
[Q, R] = qr_GSP(H)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, H)
disp(erro)

// --------------------------------------------------------
disp("Algoritmos questao 1")
disp("M1")
disp("Primeira versão")
[U, R] = qr_House_v1(M1)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M1)
disp(erro)
disp("Segunda versão")
[U, R] = qr_House_v2(M1)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M1)
disp(erro)

disp("M2")
disp("Primeira versão")
[U, R] = qr_House_v1(M2)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M2)
disp(erro)
disp("Segunda versão")
[U, R] = qr_House_v2(M2)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, M2)
disp(erro)

disp("H")
disp("Primeira versão")
[U, R] = qr_House_v1(H)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, H)
disp(erro)
disp("Segunda versão")
[U, R] = qr_House_v2(H)
[Q] = Constroi_Q_House(U)
disp("QT * Q")
disp(Q' * Q) // Tem que printar a identidade.
disp("Q * R")
disp(Q * R) // Tem que printar A.
disp("Erro: ")
erro = Compute_accuracy_QR(Q, R, H)
disp(erro)