clear
tic

co=6.8;
cm=0.5;
cu=103;

lb=[2 1 5 5 5];
ub=[20 10 50 50 50];

nl=lb(1);
fl=lb(2);
l1l=lb(3);
l2l=lb(4);
l3l=lb(5);

nu=ub(1);
fu=ub(2);
l1u=ub(3);
l2u=ub(4);
l3u=ub(5);

std0=1e100;
n0=0;
f0=0;
l10=0;
l20=0;
l30=0;


global J
global u
global p
global lam
global b
global a
global s
global R0_theta


theta=[0.124266173663911,1.54458138438586,1.01949946485847,0.484549499028519,1.70641297762755];

u=theta(1);
p=theta(2);
lam=theta(3);
b=theta(4);
a=theta(5);

J=3;
S=[45,70,95];
S0=40;%正常应力水平
SH=100;%容许上限
for k=1:J
    s(k)=(1/(S0+273.15)-1/(S(k)+273.15))/(1/(S0+273.15)-1/(SH+273.15));
end

w=30;

tm=52560;%6年

R0_u=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(sqrt(lam/u^p)*(((1-p/2)*sqrt(tm^b))+(p/2*w/u/sqrt(tm^b))));
R0_lam=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(1/lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)));
R0_b=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(lam/u^p)*log(tm)*(u*sqrt(tm^b)+w/sqrt(tm^b)));
R0_a=0;

R0_theta=[R0_u,R0_lam,R0_b,R0_a];

BB=[0,0,0,0,0,0,0;nl,fl,l1l,l2l,l3l,0,0];
bb=1;
for cb=1000:500:5000
%cb=3000;

bb=bb+1

nmax=min(nu,floor((cb-co*fl*(l1l+l2l+l3l))/(cu+(l1l+l2l+l3l)*cm)));
for n=max(nl,BB(bb,1)):nmax
    fmax=min(fu,floor((cb-n*cu-n*(l1l+l2l+l3l)*cm)/(co*(l1l+l2l+l3l))));
    for f=fl:fmax
        l1max=min(l1u,floor((cb-n*cu)/(n*cm+f*co)-l2l-l3l));
            for l1=l1l:l1max
                l2max=min(l2u,floor((cb-n*cu)/(n*cm+f*co)-l1-l3l));
                for l2=l2l:l2max
                    %l3max=min(l3u,floor((cb-n*cu)/(n*cm+f*co)-l1-l2));
                    l3=floor((cb-n*cu)/(n*cm+f*co)-l1-l2);
                    %for l3=l3l:l3max                    
                      if (l3<l3l||l3>l3u) continue              
                      end
                    std=objective_function_three([n,f,l1,l2,l3]);
                    if(std>0&&std<std0)
                    std0=std;
                    n0=n;
                    f0=f;
                    l10=l1;
                    l20=l2;
                    l30=l3;
                    end
                    end
                end
            end
    end
%end
cb0=n0*cu+n0*(l10+l20+l30)*cm+co*f0*(l10+l20+l30);
BB(bb+1,:)=[n0 f0 l10 l20 l30 std0 cb0];
end
t=toc


