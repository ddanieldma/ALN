function [x] = L_solve(L, b)
	//calcula y sendo y = Ux e Ly = b
    x=zeros(b);
    x(1)=b(1)/L(1,1)

    for i=2:n
        x(i)=(b(i)-L(i,1:(i-1))*x(1:(i-1)));
    end
endfunction