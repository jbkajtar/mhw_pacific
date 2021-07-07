function [i1,i2]=findrange(xx,x1,x2)
%
% function [i1,i2]=findrange(xx,x1,x2)
%   e.g.,: xx = lonts
%          x1 = 50
%          x2 = 70

i1=min(find(xx>=x1));
i2=max(find(xx<=x2));

