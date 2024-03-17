// 3) Modifique a função dada trocando linhas quando no início da iteração j o elemento na posição (j,j) é nulo. Chame esta nova função de Gaussian_Elimination_2 e teste-a com a matriz A1 e o vetor b1 dados. Agora teste-a com a matriz A2=[0 10-20 1; 10-20 1 1; 1 2 1] e o vetor b2=[1; 0; 0]


function [x, C]=Gaussian_Elimination_2(A, b)
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
    
    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
    end
    
    C=C(1:n,1:n);
endfunction

A2 = [0, 10^(-20), 1;
       10^(-20), 1, 1;
       1, 2, 1]
       
b2 = [1; 0; 0]

[x2, C2] = Gaussian_Elimination_2(A2, b2)

disp("Matriz C2:")
disp(C2)
disp("Vetor x2:")
disp(x2)
