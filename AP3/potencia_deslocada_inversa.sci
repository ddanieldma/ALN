function [C]=Gaussian_elimination(A, b)
    C=[A,b]
    [n]=size(C,1)
    
    for j=1:(n-1)
        // se o nosso pivo for 0
        if C(j, j) == 0 then
            // encontramos na mesma coluna, na linha a partir do pivo o primeiro elemento não nulo 
            non_zero_index = find(C(j:n, j), 1)
            // e trocamos as linhas da matriz. a expressão non_zero_index + j - 1 é porque o índice do find é relativo ao vetor que é passado para ele, não à matriz em si, então é necessário fazer um ajuste.
            C([non_zero_index + j - 1, j], :) = C([j, non_zero_index + j - 1], :)
        end
        
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1)
        end
    end
    
    C=C(1:n,1:n)
endfunction

function [x]=Resolve_com_LU(C, b)
    [n] = size(C, 1)
    L = tril(C, -1) + eye(n,n)
    U = triu(C)

    //calcula y sendo y = Ux e Ly = b
    y=zeros(b)
    y(1)=b(1)

    for i=2:n
        y(i)=(b(i)-L(i,1:(i-1))*y(1:(i-1)))
    end
    
    //calcula x sendo Ux = y
    x=zeros(y)
    x(n)=y(n)/U(n,n)
    
    for i=n-1:-1:1
        x(i)=(y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i)
    end
endfunction

function [x] = Resolve_sistema(A, b)
    // Resolve sistema linear usando a eliminação de Gauss-Seidel

    [C] = Gaussian_elimination(A, b)
    [x] = Resolve_com_LU(C, b)
endfunction

function [lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, epsilon, alfa, M)
	// Encontra o autovalor da matriz A mais próximo do alfa dado pelo método da
	// potência inversa deslocada.

	[n] = size(A, 1)
	k = 0
	x0 = x0/norm(x0, 2)
	n_erro = epsilon + 1

	I = eye(n)
	A_deslocada = A - alfa*I
	while k<=M && n_erro >= epsilon then
		[x1] = Resolve_sistema((A_deslocada), x0)
		x1 = x1/norm(x1, 2)
		lambda = x1'*A*x1
		if x1'*x0 < 0 then
			x1 = -x1
		end
		n_erro = norm((x1 - x0), 2)
		x0 = x1
		k = k+1
	end

	lambda1 = alfa + 1/max(abs(x1))
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
disp("Número de iterações: ")
disp(k)
disp("Erro das útlimas duas iterações: ")
disp(n_erro)