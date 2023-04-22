function dudx = dUdX(U)
global N dx
dudx = zeros(size(U));
for i = 2:N-1
    dudx(i) = (U(i+1) - U(i-1))/(2*dx);
end
end