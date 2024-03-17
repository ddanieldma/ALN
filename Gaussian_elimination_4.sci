// 5) Modifique a função do item 4 para escolher sempre o maior pivô em módulo no início da iteração j independente do elemento na posição (j,j) ser nulo ou não. Nessa função, retorne também a matriz de permutação P utilizada. Chame esta nova função de Gaussian_Elimination_4 e teste-a com a matriz A3 e o vetor b3 dados

// função que cria uma matriz de permutação permutando uma matriz identidade com as linhas dadas nos parâmetros
function [P] = Create_Permutation(size_m, line_1, line_2)
    //create identity matrix
    I = eye(size_m, size_m)
    // permutate its lines
    I([line_1, line_2], :) = I([line_2, line_1], :)
    P = I
endfunction

function [x, C, P]=Gaussian_Elimination_4(A, b)
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
    
    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
    end
    
    C=C(1:n,1:n);
endfunction

A3 = [10^(-20), 10^(-20), 1;
       10^(-20), 1, 1;
       1, 2, 1]
       
b3 = [1; 0; 0]

[x3, C3, P] = Gaussian_Elimination_4(A3, b3)

disp("Matriz C3:")
disp(C3)
disp("Vetor x3:")
disp(x3)
disp("Matriz de permutação:")
disp(P)
