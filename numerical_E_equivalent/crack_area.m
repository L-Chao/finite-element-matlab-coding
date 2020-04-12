function A = crack_area(R, a, b)
A = zeros(size(a));
for i = 1:length(a)
%     y0 = (a(i)*sqrt(a(i)^2*R^2 + b^4 - a(i)^2*b^2) - a(i)^2*R)/(b^2 - a(i)^2);
%     x0 = sqrt(2*R*y0 - y0^2);
    solution = fsolve(@(x) myfun(R, a(i), b, x), [R; R]);
    x0 = solution(1);
    y0 = solution(2);
    if y0 <= R
        A(i) = 2*quadgk(@(x) a(i)/b*sqrt(b^2 - x.^2) - R + sqrt(R^2 - x.^2), 0, x0);
    else
        A(i) = pi*R^2 - 2*quadgk(@(x) R + sqrt(R^2 - x.^2) - a(i)/b*sqrt(b^2 - x.^2), 0, x0);
    end
end