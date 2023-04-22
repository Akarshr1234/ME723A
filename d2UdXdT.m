function d2udxdt = d2UdXdT(Uo,Un)
global N dx dt
dudx_o = zeros(size(Uo));
dudx_n = zeros(size(Uo));
for i = 2:N-1
    dudx_o(i) = (Uo(i+1) - Uo(i-1))/(2*dx);
    dudx_n(i) = (Un(i+1) - Un(i-1))/(2*dx);
end
d2udxdt = (dudx_n - dudx_o)/dt;
end