// Utilizando a função resolve com LU e gauss elimination 4 
function [x, C, P] = Gaussian_Elimination_4(A, b)
    C = [A, b];
    [n] = size(C, 1);

    // Criamos uma matriz permutação para guardar as permutações que faremos na matriz
    P = eye(n,n);

    for j = 1:(n-1)
        linha_pra_trocar = 0;

        // Pegando o maior valor da coluna e o seu índice
        [maior_abs, idx] = max(abs(C(j:n, j)));
        // Como estamos iterando sobre j, precisamos somar o index ao j, para obter o índice correto
        idx = idx+j-1

        // Checando se linha pra trocar é diferente de 0 para assim fazer a troca
        C([j, idx], :) = C([idx, j], :);
        P([j, idx], :) = P([idx, j], :);

        for i = (j+1):n
            C(i, j) = C(i, j) / C(j, j);
            C(i, j+1:n+1) = C(i, j+1:n+1) - C(i, j) * C(j, j+1:n+1);
        end
    end

    x = zeros(n, 1);
    x(n) = C(n, n+1) / C(n, n);

    for i = n-1:-1:1
        x(i) = (C(i, n+1) - C(i, i+1:n) * x(i+1:n)) / C(i, i);
    end

    C = C(1:n, 1:n);
endfunction
    

// Função Resolve com LU
function [X] = Resolve_com_LU(C, B, P)
    // Pegando o tamanho das linhas de C
    n = size(C, 1);
    
    // Pegando o tamanho de B
    linhas_b = size(B,1)
    colunas_b = size(B,2)

    // Como nossa decomposição LU é de P*A, é importantante multiplicar P por B para evitar erro nos cálculas
    B = P*B

    // Criamos uma matriz Y de zeros
    Y = zeros(linhas_b, colunas_b);

    // Resolvemos o sistema L*Y = B
    
    // Setamos a primeira linha apenas como uma divisão, já que L é lower triangular 
    Y(1, :) = B(1, :);

    // Realizamos um loop para fazer a mesma coisa que fazemos em sistemas lineares, diminuindo dos valores que já conhecemos 
    for i = 2:n
        Y(i, :) = B(i, :) - C(i, 1:i-1) * Y(1:i-1, :);
    end

    // Agora resolveremos UX = Y, de modo análogo ao anterior, mas dessa vez devemos nos atentar ao fato de que U é upper triangular
    X = zeros(linhas_b, colunas_b);

    X(n, :) = Y(n, :) / C(n, n);
    for i = n-1:-1:1
        X(i, :) = (Y(i, :) - C(i, i+1:n) * X(i+1:n, :)) / C(i, i);
    end
endfunction


// Criando a função especificada
function [lambda1,x1,k,n_erro] = Potencia_deslocada_inversa (A,x0,epsilon,alfa,M,exibe_mensagem)    
    // Número de iterações
    k = 0;

    // Primeiros precisamos descobrir a decomposição de C da matriz (lower e upper em uma matriz só)
    // Aplicamos então a gaussian elimation 4
    tamanho_A = size(A,1)
    
    A_alpha = A - alfa * eye(tamanho_A, tamanho_A)
    [x,C, Per] = Gaussian_Elimination_4(A_alpha, zeros(size(A_alpha, 1), 1))

    // Dividimos x0 pela sua norma
    x0 = x0 / norm(x0, 2);
    // Setando o número de erro como maior que epsilon
    n_erro = epsilon + 1;
    // enquanto o número máximo de iterações não for atingido, e o erro for maior que epsilon
    while k <= M & n_erro >= epsilon then
        // Achando x1
        [x1] = Resolve_com_LU(C, x0, Per)
        
        // Transformando x1 em unitário
        x1 = x1 / norm(x1, 2)

        // Pegamos o maior elemento em módulo
        lambda1 = x1' * A * x1;

        if x1' * x0 < 0 then
            x1 = -x1;
        end

        // pegamos o erro
        n_erro = norm(x1 - x0, 2);
        // Setamos as variáveis para continuar com o loop
        x0 = x1;
        k = k + 1;
    end
    if (exibe_mensagem) then
        if k == M
            disp('Número máximo de iterações atingido. Não houve convergência');
        else
            disp('Houve Convergência para o autovalor: ', lambda1);
        end
    end
end
