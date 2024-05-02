function [P] = Create_Permutation(size_m, line_1, line_2)
    //create identity matrix
    I = eye(size_m, size_m)
    // permutate its lines
    I([line_1, line_2], :) = I([line_2, line_1], :)
    P = I
endfunction

function [C, P]=Gaussian_Elimination_4(A)
    C=[A];
    [n]=size(C,1);
    
    // matriz identidade que vai "armazenar" todas as matrizes de permutação usadas ao longo da função
    P = eye(n, n)

    for j=1:(n-1)
        // agora pegamos apenas o pivo com maior valor absoluto e usamos ele para criar a matrix de permutação e aplicá-la na nossa matriz C
        [max_value, max_index] = max(abs(C(j:n,j)));
        [Perm] = Create_Permutation(n, j, max_index + (j-1))

        // multiplicando C e a identidade pela matriz de permutação
        C = Perm * C
        P = Perm * P
        for i=(j+1):n
            C(i, j) = C(i, j) / C(j, j);
            C(i, j + 1:n) = C(i, j + 1:n) - C(i, j) * C(j, j + 1:n);
        end
    end
    
    C=C(1:n,1:n);
endfunction

function [x] = Resolve_com_LU_4(C, b, P)
    // Resolve o sistema de equações usando a matriz C e a maptriz P de
    // permutação dadas pela função G_E_4

    [n] = size(C, 1);
    L = tril(C, -1) + eye(n,n);
    U = triu(C);

    b = P*b

    //calcula y sendo y = Ux e Ly = b
    y = zeros(b);
    y(1) = b(1);

    for i = 2:n
        y(i) = (b(i) - L(i, 1:(i - 1)) * y(1:(i - 1)));
    end
    
    //calcula x sendo Ux = y
    x = zeros(y);
    x(n) = y(n) / U(n, n);
    
    for i = (n-1):-1:1
        x(i) = (y(i) - U(i, (i + 1):n) * x(i + 1:n)) / U(i, i);
    end
endfunction

function [x] = Resolve_sistema_4(A, b)
    // Resolve sistema usando a elimincação de Gauss-Seidel com trocas de
    // linha utilizando o pivô de maior módulo.

    [C, P] = Gaussian_Elimination_4(A)
    [x] = Resolve_com_LU_4(C, b, P)
endfunction

function [lambda, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, epsilon, alfa, M)
	// Encontra o autovalor da matriz A mais próximo do alfa dado pelo método da
	// potência inversa deslocada.

	// Inicializando contador.
    k = 0
    
    // Inicializando erro.
	n_erro = epsilon + 1
	
    // Normalizando o vetor inicial.
    x0 = x0/norm(x0, 2)
    
	// Criando matriz de deslocamento.
    [n] = size(A, 1)
	I = eye(n, n)
	A_deslocada = A - alfa*I

	while k<=M && n_erro >= epsilon then
        // Resolvendo sistema para encontrar próxima iteração.
        [x1] = Resolve_sistema_4((A_deslocada), x0)
		
        // Normalizando vetor.
        x1 = x1/norm(x1, 2)
		
        // Aproximando lambda.
        lambda = x1' * A * x1
		if (x1'*x0 < 0) then
			x1 = -x1
		end

        // Calculando erro da última iteração.
		n_erro = norm(x1 - x0, 2)

        // Indo para próxima iteração.
		x0 = x1
		k = k+1
	end
endfunction

// A = [0 5 -6;
// 	-4 12 -12;
// 	-2 -2 10]

// A = [1 0 1;
// 	 1 1 1;
// 	 0 1 1]
// 

A = [-2.5 3.5 -2.5;
     -1 2 1;
     -3.5 3.5 -1.5]
//  


AtA = A'*A

// x0 = [1; 1; 1]

// // Centros dos discos.
// centros = diag(AtA)

// // Número de linhas da matriz
// [n] = size(A, 1)

// for i=1:n
// 	disp("===========================================")
// 	disp("Centro do disco proximo do autovalor: ")
// 	disp(centros(i))
	
// 	[lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(AtA, x0, 0.001, centros(i), 100)
	
// 	disp("Autovalor encontrado: ")
// 	disp(lambda1)
// end

// Matriz com discos de Gerschgorin disjuntos. Como os centros estão bem espaçados
// e os raios são bem pequenos os discos são disjuntos.
A = [10 1 1/2;
     1 20 2;
     1/2 2 30]
// 

x0 = [1; 0; 0]
[n] = size(A, 1)

// Centros dos discos.
centros = diag(A)
disp("Centros dos discos: ")
disp(centros)

autovalores = zeros(n)

for i=1:n
	[lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, 0.001, centros(i), 100)
	
	autovalores(i) = lambda1
end

disp("Autovalores encontrados: ")
disp(autovalores)

disp("Autovalores dados pelo scilab: ")
disp(spec(A))