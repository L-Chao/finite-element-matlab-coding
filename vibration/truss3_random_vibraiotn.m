close all;
clear;
load('K_s.mat')
load('M_s.mat')
load('Phi.mat')
load('omega2.mat')

N = 2^11;
T = 2;
t = linspace(0, 2, N);
F = zeros(18, N);
F(1, :) = wgn(1, N, 90);
F_n = Phi'*F;
xi = 0.02;
displance = zeros(18, N);
for i = 1:length(omega2)
    displance(i, :) = responce(sqrt(omega2(i)), xi,  T, N, F(i, :));
end
d = Phi*displance;
e = (d(7:9, :) - d(16:18, :))/1.9050;
e = e(2, :);
sigma = e*6.89e10*1e-6;
%%
figure()
plot(t, F(1, :))
xlabel('t(s)','Interpreter', 'latex');
ylabel('$F(N)$', 'Interpreter', 'latex')

figure()
plot(t, sigma)
xlabel('t(s)','Interpreter', 'latex');
ylabel('$\sigma(Mpa)$', 'Interpreter', 'latex')
%%
[c,hist,edges,rmm,idx] = rainflow(sigma);
T = array2table(c,'VariableNames',{'Count','Range','Mean','Start','End'});
histogram('BinEdges',edges','BinCounts',sum(hist,2))
xlabel('Stress Range')
ylabel('Cycle Counts')