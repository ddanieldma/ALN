matriz = [1, 2, -5, 3;
          0, 5, 6, 6;
          7, 8 ,9, 9;
          4, 5, 6, 6]
          
function [P] = create_permutation(size_m, line_1, line_2)
    //create identity matrix
    I = eye(size_m, size_m)
    // permutate its lines
    I([line_1, line_2], :) = I([line_2, line_1], :)
    P = I
endfunction

[P] = create_permutation(4, 2, 3)

disp(P * matriz)
