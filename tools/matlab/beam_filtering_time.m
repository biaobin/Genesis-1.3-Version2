%% beam filtering (energy)

function [new]=beam_filtering_time(old,Q,down_time,up_time)

z = old(:,5);
np_old=length(z);
ind = find(z>=down_time*3e8 & z <= up_time*3e8);
x = old(ind,1);
px = old(ind,2);
y = old(ind,3);
py = old(ind,4);
z  = old(ind,5);
gamma = old(ind,6);
Q = Q/np_old*length(gamma);


new = struct('x',x,'px',px,'y',y,'py',py,'z',z,'gamma',gamma,'Q',Q);
end