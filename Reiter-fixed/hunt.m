function xind = hunt(xx,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hunt.m
%
% MATLAB code to bracket a linear interpolation interval.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2014)
%
% This Version : 2/26/14
% NOTE: Adapted from a similarly named routine in Numerical Recipes
% Fortran.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xind = sum(x>=xx);
    if (xind==0);
        xind=1;
    elseif (xind==length(xx));
        xind=length(xx)-1;
    end;
end