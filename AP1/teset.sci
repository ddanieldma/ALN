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

    for j=1:(n - 1)
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
    [n] = size(C, 1)
    L = tril(C, -1) + eye(n,n)
    U = triu(C)
    B = P * B
    
    disp("Perm")
    disp(B)

    //calcula y sendo y = Ux e Ly = b
    Y=zeros(B);
    Y(1,:)=B(1,:)/C(1,1);


    // for i=2:n
    //     Y(i,:)=(B(i,:)-C(i,1:i)*Y(1:i,:))/C(i,i);
    // end
    for i=2:n
        Y(i,:)=(B(i,:)-L(i,1:i)*Y(1:i,:))/L(i,i);
    end

    teste = L

    disp("teste")    
    disp(teste)
    disp(teste * Y)

    //calcula x sendo Ux = y
    X=zeros(Y);
    // Calcula x, sendo Ux=C(1:n,n+1)
    X(n,:)=Y(n,:)/C(n,n);
    // x(n)=C(n,n+1)/C(n,n);
    
    // for i=n-1:-1:1
    //     X(i,:)=(Y(i,:)-C(i,i+1:n)*X(i+1:n,:))/C(i,i);
    // end
    for i=n-1:-1:1
        X(i,:)=(Y(i,:)-U(i,i+1:n)*X(i+1:n,:))/U(i,i);
    end
endfunction

A1 = [1 1 1 0;
      1 1 0 1; 
      1 0 1 1; 
      1 1 1 1]
B1 = [2 4 -1 5;
      0 1 0 3;
      2 2 -1 1;
      0 1 1 5]
      
[C1, P1] = Gaussian_Elimination_4(A1, B1(:,1))
[X] = Resolve_com_LU(C1, B1, P1)

disp(X)
