%% beam filtering (tr_mometum)

function [new]=beam_filtering_x_y(old,Q,down_x,up_x,down_y,up_y)

x = old(:,1);
y = old(:,3);
np_old=length(x);
ind = find(x>=down_x & x <= up_x& y>=down_y & y <= up_y);
x = old(ind,1);
px = old(ind,2);
y = old(ind,3);
py = old(ind,4);
z  = old(ind,5);
gamma = old(ind,6);
Q = Q/np_old*length(gamma);


new = struct('x',x,'px',px,'y',y,'py',py,'z',z,'gamma',gamma,'Q',Q);
end