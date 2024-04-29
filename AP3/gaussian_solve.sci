function [C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
	disp("C")
	disp(C)
    [n]=size(C,1);
    
    for j=1:(n-1)
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j);
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    
    C=C(1:n,1:n);
endfunction

function [C]=Gaussian_elimination(A, b)
    C=[A,b];
    [n]=size(C,1);
    
    for j=1:(n-1)
        // se o nosso pivo for 0
        if C(j, j) == 0 then
            // encontramos na mesma coluna, na linha a partir do pivo o primeiro elemento não nulo 
            non_zero_index = find(C(j:n, j), 1);
            // e trocamos as linhas da matriz. a expressão non_zero_index + j - 1 é porque o índice do find é relativo ao vetor que é passado para ele, não à matriz em si, então é necessário fazer um ajuste.
            C([non_zero_index + j - 1, j], :) = C([j, non_zero_index + j - 1], :);
        end
        
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j);
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    
    disp(C);
    C=C(1:n,1:n);
endfunction

function [x]=Resolve_com_LU(C, b)
    [n] = size(C, 1);
    L = tril(C, -1) + eye(n,n);
    U = triu(C);

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

function [x] = Resolve_sistema(A, b)
    // Resolve sistema linear usando a eliminação de Gauss-Seidel

    [C] = Gaussian_elimination(A, b);
    [x] = Resolve_com_LU(C, b);
endfunction

A = [1 0 1;
	 0 1 1;
	 1 1 1]

A_inv = inv(A)

b = [1; 2; 3]

[x] = Resolve_sistema(A, b)

disp("x")
disp(x)

disp("Pela matriz inversa")
disp(A_inv * b)