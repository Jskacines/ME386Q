function K = K_matrix(N,E,A)

% N is number of nodes 

%% Building global stiffness matrix 

K = zeros(N,N); 

K(1,1) = 1; 
K(1,2) = -1; 
K(N,N) = 1; 
K(N,N-1) = -1; 

c = 0;

for i = 2:N-1
    for j = 2:N-1

        K(i,j-1+c) = E*A/N * -1;
        K(i,j+c) = E*A/N * 2;
        K(i,j+1+c) = E*A/N * -1;

        c = c+1;
        break
    end
end
