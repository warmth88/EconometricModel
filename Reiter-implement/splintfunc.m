function y = splintfunc(xa,ya,y2a,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% splintfunc.m
%
% MATLAB code to evaluate a function which has been approximated with 
% cubic splines using the routine splinefunc.m. Returns the interpolated
% function value at a specified input.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2014)
%
% This Version : 2/26/14
% NOTE: Adapted from a similarly named routine in Numerical Recipes
% Fortran.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    klo = hunt(xa,x);
    khi = klo+1;
    h = xa(khi)-xa(klo);
    a = (xa(khi)-x)/h;
    b = (x - xa(klo))/h;
    y = a * ya(klo) + b*ya(khi) + ((a ^ 3.0 - a)*y2a(klo)+(b^3.0-b)*y2a(khi))*(h^2.0)/6.0;


end