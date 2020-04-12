function D = differential(R, a, b, error)
D0 = 0;
H = 1e-3;
while(1)
    D = (crack_area(R, a + H, b) - crack_area(R, a - H, b))/(2*H);
    if (abs(D - D0) <= error)
        break
    end
    D0 = D;
    H = H*1/2;
end
end