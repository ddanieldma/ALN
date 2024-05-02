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

function [C]=Gaussian_elimination(A)
    C=[A];
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
            C(i,j+1:n)=C(i,j+1:n)-C(i,j)*C(j,j+1:n);
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

    [C] = Gaussian_elimination(A);
    [x] = Resolve_com_LU(C, b);
endfunction

A = [3 0 5;
	 0 5 6;
	 2 5 7]

A_inv = inv(A)

b = [1; 2; 3]

// [x] = Resolve_sistema(A, b)
[x] = Resolve_sistema_4(A, b)

disp("x")
disp(x)

disp("Pela matriz inversa")
disp(A_inv * b)