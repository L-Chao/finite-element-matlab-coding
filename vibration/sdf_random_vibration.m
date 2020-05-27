% close all;
clear;

m = 100;
k = 1e5;

omega = sqrt(k/m);
c = 632;
xi = c/sqrt(2*m*omega);

T = 0.64;
omega_bar = 2*pi/T;
M = 10;
N = 2^M;

t = linspace(0, 0.64, N);
Dt = t(2) - t(1);
tt = length(t);
F = wgn(1, N - 1, 0);
omega_n = zeros(1, N - 1);
for i = 1:tt - 1
%     F = F0(1: i);
%     if t(i) <= 0.16
%         F(i) = 12e4*t(i)/0.16;
%     elseif t(i) <= 0.48
%         F(i) = -75e4*(t(i)-0.16) + 12e4;
%     else
%         F(i) = min(0, 75e4*(t(i) - 0.64));
%     end
%     Cn = fft(F/N);
    if i <= N/2
        omega_n(i) = (i - 1)*omega_bar;
    else
        omega_n(i) = - (N - (i - 1))*omega_bar;
    end
%     rn(i) = omega_n(i)/omega;
end
rn = omega_n/omega;
Cn = fft(F/N);
uu = Cn./(k*(1 - rn.^2 + 2*xi*sqrt(-1)*rn));
u = fft(uu');
figure()
t = 0:Dt:(0.64 - Dt);
plot(t', real(u))
xlabel('Time(sec)');
ylabel('Displacement(in.)');
grid on