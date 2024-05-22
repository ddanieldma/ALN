// Quando vemos o valor de pesos, podemos ver que o primeiro peso é muito
// maior que os outros. Podemos testar a relevância, por exemplo, fazendo
// a análise de precisão do modelo zerando uma hora um peso grande e outra
// hora um peso pequeno.

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

function analise_precisao(pesos, vies, teste, teste_respostas)
    // Imprime as métricas relacionadas à precisão do modelo de preivsão.
    disp("===========================================================================")
   
    previsoes = vies + teste * pesos

    // Multiplicando cada previsão pela resposta correta.
    acertos_erros = previsoes .* teste_respostas

    acertos_erros = acertos_erros >= 0
    disp("Número de acertos")
    total = sum(acertos_erros)
    disp(total)

    // Matriz que armazena a previsão feita e se ela foi certa, para que possamos
    // fazer as contagens necessárias para a matriz de confusão.
    matriz_indicadora = [previsoes, acertos_erros]
    linhas = size(matriz_indicadora, 'r')

    // Contagem de verdadeiros positivos.
    v_p = 0
    // Contagem de verdadeiros negativos.
    v_n = 0
    // Contagem de falsos positivos.
    f_p = 0
    // Contagem de falsos negativos.
    f_n = 0

    for i = 1:linhas
        if (matriz_indicadora(i,1) >= 0 & matriz_indicadora(i,2)) then
            v_p = v_p + 1;
        elseif (matriz_indicadora(i,1) < 0 & (matriz_indicadora(i,2))) then
            v_n = v_n + 1;
        elseif (matriz_indicadora(i,1) >= 0 & ~(matriz_indicadora(i,2))) then
            f_p = f_p + 1;
        else
            f_n = f_n + 1;
        end
    end

    disp("Verdadeiros positivos:")
    disp(v_p)
    disp("Verdadeiros negativos:")
    disp(v_n)
    disp("Falsos positivos:")
    disp(f_p)
    disp("Falsos negativos:")
    disp(f_n)

    acuracia = (v_p + v_n) / total
    disp("Acurácia:")
    disp(acuracia)
    precisao = v_p / (v_p + f_p)
    disp("Precisão:")
    disp(precisao)
    recall = v_p / (v_p + f_n)
    disp("Recall")
    disp(recall)
    p_falso_alarme = f_p / (v_p + f_p)
    disp("Probabilidade de falso alarme")
    disp(p_falso_alarme)
    p_omissao_alarme = f_n / (v_n + f_n)
    disp("Probabilidade da omissão de alarme")
    disp(p_omissao_alarme)
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

// Colocando peso mais relevante positivo bem pequeno: 
pesos_original = pesos

disp("Peso mais relevante positivo pequeno")
[_, idx] = max(abs(pesos))
disp(pesos(idx))
pesos(idx) = 10**(-2);
disp("Pesos")
disp(pesos)

analise_precisao(pesos, vies, teste, teste_respostas)
// É possível perceber que ele prevê qualquer coisa como negativo.

// Colocando peso negativo mais relevante bem pequeno.
disp("Peso negativo mais relevante pequeno")
pesos = pesos_original
[_, idx] = min(pesos)
disp(pesos(idx))
pesos(idx) = 10**(-2);
disp("Pesos")
disp(pesos)

analise_precisao(pesos, vies, teste, teste_respostas)
// Falou que a grande maioria é positivo.


disp("Menor peso menor ainda")
pesos = pesos_original
[_, idx] = min(abs(pesos))
disp(pesos(idx))
pesos(idx) = 10**(-2);
disp("Pesos")
disp(pesos)

analise_precisao(pesos, vies, teste, teste_respostas)


// Colocando menor peso muito grande
disp("Menor peso muito grande")
pesos = pesos_original
[_, idx] = min(abs(pesos))
disp(pesos(idx))
pesos(idx) = 10**(2);
disp("Pesos")
disp(pesos)

analise_precisao(pesos, vies, teste,teste_respostas)


// E retirando os dois maiores pesos.
disp("Retirando dois pesos com maior relevância")
pesos = pesos_original
[_, idx] = max(abs(pesos))
disp(pesos(idx))
pesos(idx) = 10**(-2);
[_, idx] = max(abs(pesos))
disp(pesos(idx))
pesos(idx) = 10**(-2);
disp("Pesos")
disp(pesos)


analise_precisao(pesos, vies, teste,teste_respostas)