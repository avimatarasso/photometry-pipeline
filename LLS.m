function [fit405] = LLS(d470, d405)
% written by Avi, following Li's method
% Originally following Lerner TN et al. (2015) Cell. 162:635-647. 
% 

g=d470;
f=d405;

n = length(d470);
if n~=length(d405)
    warning('470 and 405 not the same length!')
end

m = (n*sum(f.*g)-sum(f)*sum(g))/(n*sum(f.^2)-(sum(f)*sum(f)));
c = (sum(g)*sum(f.^2) - (sum(f)*sum(f.*g)))/(n*sum(f.^2) - sum(f)*sum(f));

fit405 = m*d405 + c;

%if mean
end