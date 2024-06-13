function [Q, R] = qr_House(A)
	// Description of qr_House(input)

	[m, n] = size(A)
	R = A

	for j = 1:n
		x = A(j:n, j)
		x_size = size(x, 'r')
		alpha = norm(x)
		if x(1) >=0 then
			alpha = -alpha
		end

		e_1 = [1; zeros(n - j - 1, 1)]
		u = x + alpha * e_1
		u = u/norm(u)

		H_j = eye(x_size, x_size) - 2 * u * u'

		H = eye()
	end
endfunction

function [Q, R] = qr_House(A)
	
	[m, n] = size(A)
	R = A
	Q = eye(m, m)
	// k = min(m, n)
	// H = m, n

	for j = 1:n
		u = R(j:n, j)
		u_norm = norm(u)
		u_size = size(u, 'r')
		if u(1) >= 0 then
			u_norm = -u_norm
		end
		
		u(1) = u(1) - u_norm
		u = u/u_norm
		Hj = eye(u_size) - 2 * u * u'
		H = eye(m, m)
		H(j:n, j:n) = Hj
		Q = Q * H
		R = H * R
	end
	
endfunction

A = [1 1 1;
     1 0 1;
     0 1 1]

[Q, R] = qr_House(A)

disp("Q * R")
disp(Q * R)

disp("QT * Q")
disp(Q' * Q)

// import numpy as np

def householder_reflection(a):
    v = a.copy()
    v[0] += np.copysign(np.linalg.norm(a), a[0])
    v /= np.linalg.norm(v)
    H = np.eye(len(a)) - 2 * np.outer(v, v)
    return H

def qr_decomposition(A):
    n = A.shape[0]
    Q = np.eye(n)
    R = A.copy()
    
    for i in range(n - 1):
        H = np.eye(n)
        H[i:, i:] = householder_reflection(R[i:, i])
        Q = Q @ H
        R = H @ R
    
    return Q, R

# Matriz de exemplo
A = np.array([[1, 2, 4],
              [3, 8, 14],
              [2, 6, 13]], dtype=float)

Q, R = qr_decomposition(A)

// import ace_tools as tools; tools.display_dataframe_to_user(name="Matrix Q", dataframe=pd.DataFrame(Q))
// tools.display_dataframe_to_user(name="Matrix R", dataframe=pd.DataFrame(R))


