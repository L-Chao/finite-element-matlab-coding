function F = myfun(R, a, b, x)
F = [
    x(1)^2/(b^2) + x(2)^2/(a^2) - 1;
    x(1)^2 + x(2)^2 - 2*R*x(2);];
% F = [
%     a/b*sqrt(b^2 - x(1)^2) - x(2);
%     x(1) - sqrt(2*R*x(2) - x(2)^2);];
end