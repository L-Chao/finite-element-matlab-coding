function u = responce(omega, xi,  T, N, F)
omega_bar = 2*pi/T;
omega_n = zeros(1, N);
for i = 1:N - 1
    if i <= N/2
        omega_n(i) = (i - 1)*omega_bar;
    else
        omega_n(i) = - (N - (i - 1))*omega_bar;
    end
end
rn = omega_n/omega;
Cn = fft(F/N);
uu = Cn./(omega^2*(1 - rn.^2 + 2*xi*sqrt(-1)*rn));
u = real(fft(uu'));