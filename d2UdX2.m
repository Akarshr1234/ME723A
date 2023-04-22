function du2dx2 = d2UdX2(U)
global N dx
du2dx2 = zeros(size(U));
for i = 2:N-1
    du2dx2(i) = (U(i+1) - 2*U(i) + U(i-1))/(dx^2);
end
end