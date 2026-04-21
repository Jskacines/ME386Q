function F = Load_v(N)

F = zeros(N,1); 

F(1) = 1; 

F(N) = 1; 

for i = 2:N-1

    F(i) = 2;

end