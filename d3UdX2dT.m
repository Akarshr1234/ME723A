function d3udx2dt = d3UdX2dT(Uo,Un)
global N dx dt
du2dx2_o = zeros(size(Uo));
du2dx2_n = zeros(size(Uo));
for i = 2:N-1
    du2dx2_o(i) = (Uo(i+1) - 2*Uo(i) + Uo(i-1))/(dx^2);
    du2dx2_n(i) = (Un(i+1) - 2*Un(i) + Un(i-1))/(dx^2);
end
d3udx2dt = (du2dx2_n - du2dx2_o)/dt;
end
