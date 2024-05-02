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
            C(i,j)=C(i,j)/C(j,j);
            C(i,j+1:n)=C(i,j+1:n)-C(i,j)*C(j,j+1:n);
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
    y=zeros(b);
    y(1)=b(1);

    for i=2:n
        y(i)=(b(i)-L(i,1:(i-1))*y(1:(i-1)));
    end
    
    //calcula x sendo Ux = y
    x=zeros(y);
    x(n)=y(n)/U(n,n);
    
    for i=n-1:-1:1
        x(i)=(y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
    end
endfunction

function [x] = Resolve_sistema_4(A, b)
    // Resolve sistema usando a elimincação de Gauss-Seidel com trocas de
    // linha utilizando o pivô de maior módulo.

    [C, P] = Gaussian_Elimination_4(A)
    [x] = Resolve_com_LU_4(C, b, P)
endfunction

function [lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, epsilon, alfa, M)
	// Encontra o autovalor da matriz A mais próximo do alfa dado pelo método da
	// potência inversa deslocada.

	[n] = size(A, 1)
	k = 0
	x0 = x0/norm(x0, 2)
	n_erro = epsilon + 1

	I = eye(n, n)
	A_deslocada = A - alfa*I
	while k<=M && n_erro >= epsilon then
		y0 = x0
		[x1] = Resolve_sistema_4((A_deslocada), y0)
		[m, m_index] = max(abs(x1))
		if x1(m_index) < 0 then
			m = -m
		end
		y1 = x1/m
		n_erro = norm((y1 - y0), 2)
		x0 = y1
		k = k+1
	end

	lambda1 = alfa + 1/m
endfunction

A = [0 5 -6;
	-4 12 -12;
	-2 -2 10]

x0 = [1; 1; 1]

[lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, 0.001, 5, 100)

disp("Autovalor encontrado: ")
disp(lambda1)
disp("Autovetor encontrado: ")
disp(x1)
[max_value, max_index] = max(abs(x1))
if x1(max_index) < 0 then
    max_value = -max_value
end
disp(x1/max_value)
disp("Número de iterações: ")
disp(k)
disp("Erro das útlimas duas iterações: ")
disp(n_erro)