cd("C:\Users\danie\OneDrive\Área de Trabalho\ALN\sailebi")

// Reading data.
cobb_douglas_data = read_csv('.\cobb_douglas.csv')

// Organizing the least squares.
// Getting the last two columns of the csv file.
A = cobb_douglas_data(:,2:3)
ones_arr = ones((A(:,1)))

disp(A)
disp(ones_arr)