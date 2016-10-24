function y2 = splinefunc(x,y,yp1,ypn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% splinefunc.m
%
% MATLAB code to approximate a function with cubic splines.  Uses a 
% natural spline endpoint condition. Returns a second-derivative array
% which can be used with splintfunc.m to evaluate the spline approx.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2014)
%
% This Version : 2/26/14
% NOTE: Adapted from a similarly named routine in Numerical Recipes
% Fortran.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(x);
u=zeros(n,1);
y2=zeros(n,1);

if (yp1>0.99e30);
    y2(1) = 0.0;
    u(1) = 0.0;
else
    y2(1) = -0.5;
    u(1) = (3.0/(x(2) - x(1))) * ( (y(2) - y(1)) / (x(2) - x(1)) - yp1);
end;

for i=2:(n-1);
        sig = (x(i) - x(i-1))/(x(i+1)-x(i-1));
        p = sig * y2(i-1) + 2.0;
        y2(i)=(sig-1.0)/p;
        u(i) = (6.0 * ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))...
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig * u(i-1))/p;
end

if (ypn>0.99e30);
    qn=0.0;
    un=0.0;
else
    qn = 0.5;
    un = (3.0/(x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)));
end;

y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0);
for k=(n-1:-1:1);
    y2(k)=y2(k)*y2(k+1)+u(k);
end

end