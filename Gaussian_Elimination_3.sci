// 4) Modifique a função do item 3 para escolher o maior pivô em módulo quando no início da iteração j o elemento na posição (j,j) é nulo. Chame esta nova função de Gaussian_Elimination_3 e teste-a com a matriz A2 e o vetor b2 dados. Agora com a matriz A3=[10-20 10-20 1; 10-20 1 1; 1 2 1] e o vetor b3=b2

function [x, C]=Gaussian_Elimination_3(A, b)
    C=[A,b];
    [n]=size(C,1);
    
    for j=1:(n-1)
        // se o nosso pivo for 0
        if C(j, j) == 0 then
            // encontramos na mesma coluna, na linha a partir do pivo o maior elemento, salvando o seu índice 
            [max_value, max_index] = max(abs(C(j:n,:)));
            C([max_index + j - 1, j], :) = C([j, max_index + j - 1], :);
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

mat = [10^(-20), 10^(-20), 1;
       10^(-20), 1, 1;
       1, 2, 1]
       
vec = [1; 0; 0]

[x, C] = Gaussian_Elimination_3(mat, vec)

disp(C)
disp(x)

disp(mat * x)
