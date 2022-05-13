%% beam filtering (energy)

function [new]=beam_filtering_gamma(old,Q,down_gamma,up_gamma)

gamma = old(:,6);
np_old=length(gamma);
ind = find(gamma>=down_gamma & gamma <= up_gamma);
x = old(ind,1);
px = old(ind,2);
y = old(ind,3);
py = old(ind,4);
z  = old(ind,5);
gamma = old(ind,6);
Q = Q/np_old*length(gamma);


new = struct('x',x,'px',px,'y',y,'py',py,'z',z,'gamma',gamma,'Q',Q);
end