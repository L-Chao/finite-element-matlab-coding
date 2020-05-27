close all;
clear;

%%%%-given values-
m = 100;
k = 1e5;

omega = sqrt(k/m);%natural frequency
xi = 0.1;%damping ratio
c = 2*m*omega*xi;%damping coefficient

T = 0.64;
omega_bar = 2*pi/T;%Excitation frequency;
M = 3;
N = 2^M;%the number of time increments N

t = 0:0.08:0.64;
Dt = 0.08;
tt = length(t);

for i = 1:tt - 1
    if t(i) <= 0.16
        F(i) = 12e4*t(i)/0.16;
    elseif t(i) <= 0.48
        F(i) = -75e4*(t(i) - 0.16) + 12e4;
    else
        F(i) = min(0, 75e4*(t(i) - 0.64));
    end
    Cn = fft(F/N);
    if i < N/2
        omega_n(i) = (i - 1)*omega_bar;
    else
        omega_n(i) = -(N - (i - 1))*omega_bar;
    end
    rn(i) = omega_n(i)/omega;
end
uu = Cn./(k*(1 - rn.^2 + 2*xi.*sqrt(-1)*rn));
u = fft(uu');
figure()
t = 0:Dt:(0.64 - Dt);
plot(t', real(u))
xlabel('Time(sec)');
ylabel('Displacement(in.)');
grid on