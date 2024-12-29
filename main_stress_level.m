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
a=theta(5);




w=30;

tm=52560;%6年

R0=1-normcdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1);
R0_u=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(sqrt(lam/u^p)*(((1-p/2)*sqrt(tm^b))+(p/2*w/u/sqrt(tm^b))));
R0_lam=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(1/lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)));
R0_b=-normpdf(sqrt(lam/u^p)*(u*sqrt(tm^b)-w/sqrt(tm^b)),0,1)*(1/2*sqrt(lam/u^p)*log(tm)*(u*sqrt(tm^b)+w/sqrt(tm^b)));
R0_a=0;
R0_theta=[R0_u,R0_lam,R0_b,R0_a];


BB=[NaN,NaN,NaN,NaN,NaN,NaN;nl,fl,l1l,l2l,NaN,NaN];
XX=NaN;
YY=NaN;
bb=1;
%for cb=1000:500:5000
cb=1000;


%计算Fisher信息矩阵
%先计算两个应力水平的情形
J=2;

% S1=40:1:99;
% S2=41:1:100;

S1=40:10:90;
S2=50:10:100;

for i=1:length(S1)
    for j=1:length(S2)
        if S1(i)>=S2(j) continue
        end

bb=bb+1;

S=[S1(i),S2(j)];
S0=40;%正常应力水平
SH=100;%容许上限
for k=1:J
    s(k)=(1/(S0+273.15)-1/(S(k)+273.15))/(1/(S0+273.15)-1/(SH+273.15));
end


std0=1e100;
n0=0;
f0=0;
l10=0;
l20=0;

nmax=min(nu,floor((cb-co*fl*(l1l+l2l))/(cu+(l1l+l2l)*cm)));
for n=nl:nmax
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

XX(bb+1)=S1(i);
YY(bb+1)=S2(j);
%end
    end
end



subplot(1,2,1)
%%%绘制渐近方差关于S1 S2的三维图
[x,y] = meshgrid(S1,S2);
bbb=2;
z=NaN(length(S1),length(S2));
for i=1:length(S1)
    for j=1:length(S2)
        if S1(i)>=S2(j) continue
        end
        bbb=bbb+1;
        if bbb>length(XX) break
        end
        z(i,j)=BB(bbb,5);
    end
end


mesh(x(1:1:length(S1),1:1:length(S1)),y(1:1:length(S1),1:1:length(S1)),z(1:1:length(S1),1:1:length(S1))')
%axis square; %使三个坐标轴等长
%grid on%添加网格 off丢掉网格线
colormap jet
colorbar
%shading interp 
xlabel('Stress level S_1')
ylabel('Stress level S_2')
zlabel('Standard deviation')




subplot(1,2,2)
%%固定S1=40:1:99，绘制标准差关于S2-S1的变化规律
for kk=1:length(S1)
    Xnew3=YY(XX==S1(kk))-XX(XX==S1(kk));
    Ynew3=BB(XX==S1(kk),5);
    if S1(kk)>=40&S1(kk)<50%if S1(kk)==40|S1(kk)==45%
       p1=plot(Xnew3,Ynew3,'r*-','LineWidth',0.5,'markersize',5)
    end
    hold on
    if S1(kk)>=50&S1(kk)<60%if S1(kk)==50|S1(kk)==55%
       p2=plot(Xnew3,Ynew3,'g+-','LineWidth',0.5,'markersize',5)
    end
    if S1(kk)>=60&S1(kk)<70%if S1(kk)==60|S1(kk)==65%
       p3=plot(Xnew3,Ynew3,'bo-','LineWidth',0.5,'markersize',5)
    end
    if S1(kk)>=70&S1(kk)<80%if S1(kk)==70|S1(kk)==75%
       p4=plot(Xnew3,Ynew3,'csquare-','LineWidth',0.5,'markersize',5)
    end
    if S1(kk)>=80&S1(kk)<90%if S1(kk)==80|S1(kk)==85%
       p5=plot(Xnew3,Ynew3,'m^-','LineWidth',0.5,'markersize',5)
    end
    if S1(kk)>=90&S1(kk)<100%if S1(kk)==90|S1(kk)==95%
       p6=plot(Xnew3,Ynew3,'kpentagram-','LineWidth',0.5,'markersize',5)
    end
end
xlabel('S_2-S_1')
ylabel('Standard deviation')
%legend([p1 p2 p3 p4 p5 p6],{'40\leq S_1<50','50\leq S_1<60','60\leq S_1<70','70\leq S_1<80','80\leq S_1<90','90\leq S_1<100'})
legend([p1 p2 p3 p4 p5 p6],{'S_1=40:1:49','S_1=50:1:59','S_1=60:1:69','S_1=70:1:79','S_1=80:1:89','S_1=90:1:99'})
%legend('Location','Northwest')

t=toc