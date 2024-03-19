function [P] = Create_Permutation(size_m, line_1, line_2)
    //create identity matrix
    I = eye(size_m, size_m)
    // permutate its lines
    I([line_1, line_2], :) = I([line_2, line_1], :)
    P = I
endfunction

function [C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    [n]=size(C,1);
    
    // matriz identidade que vai "armazenar" todas as matrizes de permutação usadas ao longo da função
    P = eye(n, n)

    for j=1:(n-1)
        // agora pegamos apenas o pivo com maior valor absoluto e usamos ele para criar a matrix de permutação e aplicá-la na nossa matriz C
        [max_value, max_index] = max(abs(C(j:n,:)));
        [Perm] = Create_Permutation(n, j, max_index)

        // multiplicando C e a identidade pela matriz de permutação
        C = Perm * C
        P = Perm * P
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j);
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    
    C=C(1:n,1:n);
endfunction

function [X]=Resolve_com_LU(C, B, P)
    [n] = size(B, 1)
    [m] = size(B, 2)
    L = tril(C, -1) + eye(n,n)
    U = triu(C)
    B = P * B

    //calcula y sendo y = Ux e Ly = b
    Y=zeros(B);
    Y(1,:)=B(1,:)

    for i=2:n
        Y(i,:)=(B(i,:)-L(i,1:(i-1))*Y(1:(i-1),:));
    end

    teste = Y
    
    disp("teste")
    disp("matriz B permutada")
    disp(B)
    disp("matriz Y")
    disp(Y)
    disp("L * Y, que deveria ser igual a B")
    disp(L * Y)
    
    //calcula x sendo Ux = y
    X=zeros(Y);
    X(n,:)=Y(n,:)/U(n,n);
    
    for i=n-1:-1:1
        X(i,:)=(Y(i,:)-U(i,i+1:m)*X(i+1:n,:))/U(i,i);
    end
endfunction

A1 = [1 -2 5 0;
      2 -4 1 3; 
      -1 1 0 2; 
      0 3 3 1]
B1 = [2 4 -1 5;
      0 1 0 3;
      2 2 -1 1;
      0 1 1 5]
      
[C1, P1] = Gaussian_Elimination_4(A1, B1(:,1))
[X] = Resolve_com_LU(C1, B1, P1)

disp(X)

A2 = [0, 10^(-20), 1;
       10^(-20), 1, 1;
       1, 2, 1]
       
B2 = [1 1 2;
      1 -1 0;
      1 0 1]
      
[C2, P2] = Gaussian_Elimination_4(A2, B2(:,1))
[X] = Resolve_com_LU(C2, B2, P2)

disp("X e A2 * X")
disp(X)
disp(A2 * X)