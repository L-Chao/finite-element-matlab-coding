clear
[X, Y] = meshgrid(0:0.2:1, [0.1; 0.2; 0.4; 0.6; 0.8; 1.0;]);
V = [
    1.0501, 1.0698, 1.1803, 1.4217, 1.8490, 2.5990;
    0.9859, 1.0123, 1.1225, 1.3542, 1.7645, 2.4871;
    0.9126, 0.9423, 1.0460, 1.2575, 1.6318, 2.2918;
    0.8303, 0.8599, 0.9506, 1.1314, 1.4510, 2.0133;
    0.7388, 0.7651, 0.8365, 0.9761, 1.2221, 1.6514;
    0.6383, 0.6578, 0.7036, 0.7915, 0.9450, 1.2062;];
ab_mat = [
    2, 5;
    5, 10;
    10, 20;
    15, 30;
    17.5, 35;
    20, 40;
    25, 40;
    30, 40;]*1e-3;
R = 25e-3;
L = 426e-3;
nu = 0.33;
A_section = pi*R^2;
E_ratio = zeros(size(ab_mat, 1), 1);
A_0 = zeros(size(ab_mat, 1), 1);
A_a = zeros(size(ab_mat, 1), 1);
for i = 1:size(ab_mat, 1)
    a_c = ab_mat(i, 1);
    b = ab_mat(i, 2);
    A_0(i) = crack_area(R, a_c, b);
    A_a(i) = differential(R, a_c, b, 1e-8);
    if i == 1
        f = 0.9126;
    elseif i == 8
        f = 1.6514;
    else
        f = interp2(X, Y, V, a_c/b, a_c/R);
    end
    integrand = @(a) (f*sqrt(pi*a)).^2.*differential(R, a, b, 1e-8);
    E_ratio(i) = 1/((1 - A_0(i)/A_section)^(-1) + (L*A_section/(2*(1 - nu^2)*quadgk(integrand, 0, a_c)))^(-1));
end
