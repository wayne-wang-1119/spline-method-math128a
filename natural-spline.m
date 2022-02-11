function [b, c, d] = ncspline(x, a)
x = x(:)'; a = a(:)'; % Make sure inputs are row vectors
h = diff(x);
n = length(x)-1; % Form tridiagonal matrix for clamped cubic splines A = sparse(2:n, 1:n-1, h(1:n-1), n+1, n+1) + ...
sparse(2:n, 3:n+1, h(2:n), n+1, n+1) + ...
sparse(2:n, 2:n, 2*(h(1:n-1) + h(2:n)), n+1, n+1); A(1,1) = 1; A(n+1,n+1) = 1;
% Note: This is the b in the right-hand side of Ax = b in Burden & Faires b = [0, 3./h(2:n).*(a(3:n+1) - a(2:n)) - 3./h(1:n-1).*(a(2:n) - a(1:n-1)), 0]'; c = (A\b)'; % Note: This is the b in the actual spline coefficients
b = (a(2:n+1) - a(1:n))./h - h./3.*(2*c(1:n) + c(2:n+1));
d = (c(2:n+1) - c(1:n))./(3*h); c = c(1:n);
end
function yy = splineeval(x, a, b, c, d, xx)
n = length(x)-1; % Number of spline segments yy = 0*xx; % Output vector
for i = 1:n % Find points in this segment
ix = xx >= x(i) & xx <= x(i+1); % Evaluate polynomial yy(ix) = a(i) + ...
b(i)*(xx(ix)-x(i)) + ...
c(i)*(xx(ix)-x(i)).^2 + ...
d(i)*(xx(ix)-x(i)).^3;
end end
t = [0,1,2,3,4,5];
a = [1.0,1.5,2.0,2.0,2.5,2.5]; % Fit natural cubic spline [b1,c1,d1] = ncspline(t, a);
xx = linspace(t(1), t(end), 1000);
yy1 = splineeval(t, a, b1, c1, d1, xx)
t = [0,1,2,3,4,5];
a1 = [1.0,0.5,1.0,1.5,1.5,1.0]; % Fit natural cubic spline [b2,c2,d2] = ncspline(t, a1);
yy2 = splineeval(t, a1, b2, c2, d2, xx)
