function [x, C]=Gaussian_Elimination_1(A, b)
    C=[A,b];
    [n]=size(C,1);
    
    for j=1:(n-1)
     //O pivô está na posição (j,j)
        for i=(j+1):n
          //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i  Linha i - C(i,j)*Linha j
            //Somente os elementos da diagonal ou acima dela são computados
            //(aqueles que compõem a matrix U)
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


// 1) Teste a função dada usando algumas matrizes quadradas A e respectivos vetores b. Use exemplos dos quais você saiba a resposta para verificar se a função realmente está funcionando corretamente.
disp("Primeiro exeemplo")
mat_1 = [1,1,0;
         0,1,1;
         1,0,1];
vec_1 = [1;1;1];

disp(mat_1)
disp(vec_1)

[x_1, lu_1] = Gaussian_Elimination_1(mat_1, vec_1);

disp(lu_1)
disp(x_1)
// É possível perceber que as matrizes L e U ficam "conjugadas" em uma única matriz, que, se lermos da diagonal para cima ou da diagonal para baixo, temos as matrizes L e U separadamente.

disp("========================")

disp("Segundo exemplo")
mat_1 = [1,2,3;
         4,-2,6;
         1,4,5];
vec_1 = [1;1;1];

disp(mat_1)
disp(vec_1)

[x_1, lu_1] = Gaussian_Elimination_1(mat_1, vec_1);

disp(lu_1)
disp(x_1)


// 2) Agora teste com a matriz A1=[1 -2 5 0; 2 -4 1 3; -1 1 0 2; 0 3 3 1] e com o vetor b1=[1;0;0;0]
matriz_q2 = [1 -2 5 0;
             2 -4 1 3; 
             -1 1 0 2; 
             0 3 3 1]
vetor_q2 = [1; 0; 0; 0]
