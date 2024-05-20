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


cd("C:\Users\danie\OneDrive\Área de Trabalho\ALN\sailebi")

// Lendo dados.
cobb_douglas_data = csvRead('cobb_douglas.csv', ";")

disp(cobb_douglas_data)

// Organizando os mínimos quadrados.
// Pegando as últimas duas colunas dos dados.
A = cobb_douglas_data(:,2:3)
ones_arr = ones((A(:,1)))

// Aplicando log antes de adicionar os 1's.
A = log(A)
A = [ones_arr, A]

disp(ones_arr)

// Criando vetor b com os "valores esperados".
b = cobb_douglas_data(:, 1)

// Logarítmo de A e b.
b = log(b)

disp(A)
disp(b)

At = A'
AtA = At * A

disp(AtA)
disp(At * b)

[x_barra] = Resolve_sistema_4(AtA, At * b)

disp(x_barra)

b_eq = x_barra(1)
b_eq = exp(b_eq)
alfa = x_barra(2)

disp("b:")
disp(b)
disp("alfa:")
disp(alfa)

// Valor pra 1910
disp("Valor para 1910 (deve ser 159)")
prev_1910 = b_eq * 147**alfa * 208**(1 - alfa)
disp(prev_1910)
disp("Erro: ")
disp(abs(159 - prev_1910))

// Comentário: conseguimos um erro de 10 em relação ao número de fato.

// Valor pra 1920
disp("Valor para 1920 (deve ser 231)")
prev_1920 = b_eq * 194**alfa * 407**(1 - alfa)
disp("Erro: ")
disp(abs(231 - prev_1920))