function [a_, b_, miu_, sigma_, e_, f_, L_]=V_I(T,X,y)
%Initialize the parameters: a0,b0,e0,f0
sizex=size(X);
N=sizex(1);
d=sizex(2);
a0=10^(-16);a_=a0;
b0=10^(-16);
b_=zeros(1,d);
for k=1:d
    b_(1,k)=b0;
end
e_=1;e0=1;
f_=1;f0=1;

%Calculate some necessary terms in regard to X and y
XXT=zeros(d,d);
yX=zeros(d,1);
yy=0;
yXT=zeros(1,d);
XTX=0;
for i=1:N
    XXT=XXT+X(i,:)'*X(i,:);
    yX=yX+X(i,:)'*y(i,1);
    yy=yy+y(i,1)^2;
    yXT=yXT+y(i,1)*X(i,:);
    XTX=XTX+X(i,:)*X(i,:)';    
end
L_=zeros(1,T);



%The VI iteration
for t=1:T
    %Update parameter of q(w): miu_ (d*1 vector), sigma_ (d*d matrix)
    ab_diag=zeros(d,d);
    for i=1:d
        for k=1:d
            if i==k
                ab_diag(i,k)=a_/b_(1,k);
            end
        end
    end
    miu_=(e_/f_)*((e_*XXT/f_+ab_diag)^(-1))*yX;
    sigma_=(e_*XXT/f_+ab_diag)^(-1);
    %Update parameter of q(lambda): e_, f_
    e_=N/2+e0;
    sumw=0;
    for i=1:N 
        sumw=sumw+(y(i,1)-X(i,:)*miu_)^2+X(i,:)*sigma_*X(i,:)';
    end
    f_=f0+0.5*sumw;
    %Update parameter of q(alpha):a_, b_=[b_1,...b_d]
    miumiuT=miu_*miu_';
    a_=a0+0.5;
    for k=1:d
    b_(1,k)=b0+0.5*(miumiuT(k,k)+sigma_(k,k));
    end
    
    %Calculate L function:
    %Calculate Eq(lnP)
    lnb=0;
    sum=0;
    for k=1:d
        lnb=lnb+log(b_(1,k));
        sum=sum+((miumiuT(k,k)+sigma_(k,k))/2-b0)*(a_/b_(1,k));
    end
    E=-(N+d)*log(2*pi)/2+d*a0*log(b0)-d*gammaln(a0)+e0*log(f0)-gammaln(e0)+(-yy/2+yXT*miu_-0.5*trace(XXT*sigma_)-0.5*miu_'*XXT*miu_-f0)*(e_/f_)+(a0-0.5)*(d*psi(a_)-lnb)-sum+(e0-1+N/2)*(psi(e_)-log(f_));
    %Calculate Eq(w)[lnq(w)]
    W=-d*log(2*pi)/2-0.5*log_det(sigma_)-d/2;
    %Calculate Eq(lambda)[lnq(lambda)]
    L=-gammaln(e_)+(e_-1)*psi(e_)+log(f_)-e_;
    A=-d*gammaln(a_)+d*(a_-1)*psi(a_)-d*a_+lnb;
    L_(1,t)=E-W-L-A;
end
end