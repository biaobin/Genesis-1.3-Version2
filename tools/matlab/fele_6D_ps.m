
function [phase]=fele_6D_ps(x1)





c=3*10^8;







A=textread(x1, '' , 'headerlines', 2,'commentstyle', 'shell');
A(:,5)= -A(:,5);
A(:,5)=(A(:,5)-min(A(:,5)))*c;
phase=[A(:,1),A(:,2),A(:,3),A(:,4),A(:,5),A(:,6)];



%% emit
x=A(:,1);
x_=A(:,2);
x=x-mean(x);
x_=x_-mean(x_);
emitx=sqrt(mean(x.^2)*mean(x_.^2)-(mean(x.*x_))^2)
clear x x_;
y=A(:,3);
y_=A(:,4);
y=y-mean(y);
y_=y_-mean(y_);
emity=sqrt(mean(y.^2)*mean(y_.^2)-(mean(y.*y_))^2)
clear y y_;




%delete(x1);

end