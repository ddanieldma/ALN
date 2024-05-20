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
cancer_train_data = csvRead('cancer_train_2024.csv', ";")
cancer_test_data = csvRead('cancer_test_2024.csv', ";")

// Encontrando hiperplano que faz as previsões.
// Removendo última coluna para fazer a matriz A.
colunas_a = size(cancer_train_data, 'c')
A = cancer_train_data(:,1:(colunas_a - 1))
// Vetor b é a última coluna do csv.
b = cancer_train_data(:,colunas_a)

// Adicionando coluna de 1's em A para os vieses.
A_ones = ones(A(:,1))
A = [A_ones, A]

At = A'
AtA = At * A

// Vetor com os parâmetros do hiperplano.
x_barra = Resolve_sistema_4(AtA, At * b)

// Separando x_barra entre vies e pesos.
vies = x_barra(1,:)
pesos = x_barra(2:size(x_barra, 'r'),:)

// Prevendo com todo o conjunto de dados.
// Removendo a última coluna, que contém as respostas.
colunas_teste = size(cancer_test_data, 'c')
teste = cancer_test_data(:,1:(colunas_teste - 1))
teste_respostas = cancer_test_data(:,colunas_teste)

previsoes = vies + teste * pesos

// Multiplicando cada previsão pela resposta correta.
acertos_erros = previsoes .* teste_respostas

acertos_erros = acertos_erros >= 0
disp("Número de acertos")
disp(sum(acertos_erros))

// Matriz que armazena a previsão feita e se ela foi certa, para que possamos
// fazer as contagens necessárias para a matriz de confusão.
matriz_indicadora = [previsoes, acertos_erros]
linhas = size(matriz_indicadora, 'r')

// Contagem de verdadeiros positivos.
true_p = 0
// Contagem de verdadeiros negativos.
true_n = 0
// Contagem de falsos positivos.
false_p = 0
// Contagem de falsos negativos.
false_n = 0


for i = 1:linhas
	if (matriz_indicadora(i,1) >= 0 & matriz_indicadora(i,2)) then
		true_p = true_p + 1;
	elseif (matriz_indicadora(i,1) < 0 & (matriz_indicadora(i,2))) then
		true_n = true_n + 1;
	elseif (matriz_indicadora(i,1) >= 0 & matriz_indicadora(i,2)) then
		false_p = false_p + 1;
	else
		false_n = false_n + 1;
	end
end

disp("Verdadeiros positivos:")
disp(true_p)
disp("Verdadeiros negativos:")
disp(true_n)
disp("Falsos positivos:")
disp(false_p)
disp("Falsos negativos:")
disp(false_n)