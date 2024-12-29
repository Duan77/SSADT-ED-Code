clear
tic

global J
global u
global p
global lam
global b
global a
global s
global R0_theta

co=6.8;
cm=0.5;
cu=103;

lb=[2 1 5 5];
ub=[20 10 50 50];

nl=lb(1);
fl=lb(2);
l1l=lb(3);
l2l=lb(4);

nu=ub(1);
fu=ub(2);
l1u=ub(3);
l2u=ub(4);

theta=[0.124266173663911,1.54458138438586,1.01949946485847,0.484549499028519,1.70641297762755];

u=theta(1);
p=theta(2);
lam=theta(3);
b=theta(4);
a=theta(5);%*(1+2.5/100)


%计算Fisher信息矩阵
%先计算两个应力水平的情形
J=2;
S=[45,95];
S0=40;%正常应力水平
SH=100;%容许上限
for k=1:J
    s(k)=(1/(S0+273.15)-1/(S(k)+273.15))/(1/(S0+273.15)-1/(SH+273.15));
end

w=30;

tm=52560;%6年

R0=1-normcdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1);
R0_u=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(sqrt(lam/u^p)*(((1-p/2)*sqrt(tm^b))+(p/2*w/u/sqrt(tm^b))));
R0_lam=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(1/lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)));
R0_b=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(lam/u^p)*log(tm)*(u*sqrt(tm^b)+w/sqrt(tm^b)));
R0_a=0;
R0_theta=[R0_u,R0_lam,R0_b,R0_a];

BB=[0,0,0,0,0,0;nl,fl,l1l,l2l,0,0];
bb=1;
%for cb=1000:500:5000
cb=1000;

bb=bb+1

std0=1e100;
n0=0;
f0=0;
l10=0;
l20=0;

nmax=min(nu,floor((cb-co*fl*(l1l+l2l))/(cu+(l1l+l2l)*cm)));
for n=max(nl,BB(bb,1)):nmax
    fmax=min(fu,floor((cb-n*cu-n*(l1l+l2l)*cm)/(co*(l1l+l2l))));
    for f=fl:fmax
        l1max=min(l1u,floor((cb-n*cu)/(n*cm+f*co)-l2l));
            for l1=l1l:l1max
                %l2max=min(l2u,floor((cb-n*cu)/(n*cm+f*co)-l1));
                l2=floor((cb-n*cu)/(n*cm+f*co)-l1);
                %for l2=l2l:l2max
                 if (l2<l2l||l2>l2u) continue
                 %TCost=n*cu+n*cm*(l1+l2)+f*co*(l1+l2);
                 %if (TCost<0.99*cb||TCost>cb) continue
                end
                std=objective_function_two([n,f,l1,l2]);
                if(std>0&&std<std0)
                    std0=std;
                    n0=n;
                    f0=f;
                    l10=l1;
                    l20=l2;
                end
                %end%
            end
        end
    end
%end
cb0=n0*cu+n0*(l10+l20)*cm+co*f0*(l10+l20);
BB(bb+1,:)=[n0 f0 l10 l20 std0 cb0];
%end
t=toc


