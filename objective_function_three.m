function Std=objective_function_three(x)

n=x(1);
f=x(2);
l1=x(3);
l2=x(4);
l3=x(5);

global J
global u
global p
global lam
global b
global a
global s
global R0_theta

%I=zeros(5);%参数分别为mu,p,lambda,beta,alpha
tau=[0,0,f*l1,f*(l1+l2)];
%t=[0:f:(f*l1),((f*l1)+f):f:(f*l1+f*l2),((f*l1+f*l2)+f):f:(f*l1+f*l2+f*l3)];
t=[0:f:(f*l1+f*l2+f*l3)];
ll=[0,l1,l1+l2,l1+l2+l3];
ss=[0,0,s(1),s(2),s(3)];

for k=1:J
    P1((1+ll(k)):ll(k+1))=(t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))+exp(a*ss(k)/b)*tau(k);
    P2((1+ll(k)):ll(k+1))=(t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))+exp(a*ss(k)/b)*tau(k);
    P1_b((1+ll(k)):ll(k+1))=-((t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)*a*ss(k+2)/b^2+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*a*ss(k+1)/b^2+exp(a*ss(k)/b)*tau(k)*a*ss(k)/b^2);
    P2_b((1+ll(k)):ll(k+1))=-((t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)*a*ss(k+2)/b^2+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*a*ss(k+1)/b^2+exp(a*ss(k)/b)*tau(k)*a*ss(k)/b^2);
    P1_a((1+ll(k)):ll(k+1))=(t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)*ss(k+2)/b+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*ss(k+1)/b+exp(a*ss(k)/b)*tau(k)*ss(k)/b;
    P2_a((1+ll(k)):ll(k+1))=(t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)*ss(k+2)/b+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*ss(k+1)/b+exp(a*ss(k)/b)*tau(k)*ss(k)/b;
    P1_b2((1+ll(k)):ll(k+1))=(t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)*a*ss(k+2)*(a*ss(k+2)/b^4+2/b^3)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*a*ss(k+1)*(a*ss(k+1)/b^4+2/b^3)+exp(a*ss(k)/b)*tau(k)*a*ss(k)*(a*ss(k)/b^4+2/b^3);
    P2_b2((1+ll(k)):ll(k+1))=(t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)*a*ss(k+2)*(a*ss(k+2)/b^4+2/b^3)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*a*ss(k+1)*(a*ss(k+1)/b^4+2/b^3)+exp(a*ss(k)/b)*tau(k)*a*ss(k)*(a*ss(k)/b^4+2/b^3);
    P1_b_a((1+ll(k)):ll(k+1))=-((t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)*ss(k+2)/b^2*(a*ss(k+2)/b+1)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*ss(k+1)/b^2*(a*ss(k+1)/b+1)+exp(a*ss(k)/b)*tau(k)*ss(k)/b^2*(a*ss(k)/b+1));
    P2_b_a((1+ll(k)):ll(k+1))=-((t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)*ss(k+2)/b^2*(a*ss(k+2)/b+1)+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*ss(k+1)/b^2*(a*ss(k+1)/b+1)+exp(a*ss(k)/b)*tau(k)*ss(k)/b^2*(a*ss(k)/b+1));
    P1_a2((1+ll(k)):ll(k+1))=(t((2+ll(k)):(ll(k+1)+1))-tau(k+1))*exp(a*ss(k+2)/b)*(ss(k+2)/b)^2+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*(ss(k+1)/b)^2+exp(a*ss(k)/b)*tau(k)*(ss(k)/b)^2;
    P2_a2((1+ll(k)):ll(k+1))=(t((1+ll(k)):(ll(k+1)))-tau(k+1))*exp(a*ss(k+2)/b)*(ss(k+2)/b)^2+exp(a*ss(k+1)/b)*(tau(k+1)-tau(k))*(ss(k+1)/b)^2+exp(a*ss(k)/b)*tau(k)*(ss(k)/b)^2;
    rou((1+ll(k)):ll(k+1))= P1((1+ll(k)):ll(k+1)).^b- P2((1+ll(k)):ll(k+1)).^b;
    for j=((1+ll(k)):ll(k+1))
        E_pow(j)=integral(@(y)beijihanshu_pow2p(y,rou(j),u,p,lam),0,Inf,'RelTol',1e-16,'AbsTol',1e-16);
    end
end
E_d=2*((E_pow-u^(2-p))/((1-p)*(2-p)));
E_d_rou=2*((u^(2-p)-E_pow)./((1-p)*rou));
E_d_rou2=2*((2*u^(2-p)+(p-3)*E_pow)./((p-1)*rou.^2));

rou_b(1)=P1(1).^b.*(log(P1(1))+b./P1(1).*P1_b(1));
rou_b(2:(l1+l2+l3))=P1(2:(l1+l2+l3)).^b.*(log(P1(2:(l1+l2+l3)))+b./P1(2:(l1+l2+l3)).*P1_b(2:(l1+l2+l3)))-P2(2:(l1+l2+l3)).^b.*(log(P2(2:(l1+l2+l3)))+b./P2(2:(l1+l2+l3)).*P2_b(2:(l1+l2+l3)));
rou_a(1)=b*P1(1).^(b-1).*P1_a(1);
rou_a(2:(l1+l2+l3))=b*P1(2:(l1+l2+l3)).^(b-1).*P1_a(2:(l1+l2+l3))-b*P2(2:(l1+l2+l3)).^(b-1).*P2_a(2:(l1+l2+l3));
rou_b2(1)=P1(1).^b*((log(P1(1))+b./P1(1)*P1_b(1))^2+2./P1(1).*P1_b(1)-b./P1(1).^2.*P1_b(1).^2+b./P1(1).*P1_b2(1));
rou_b2(2:(l1+l2+l3))=P1(2:(l1+l2+l3)).^b.*((log(P1(2:(l1+l2+l3)))+b./P1(2:(l1+l2+l3)).*P1_b(2:(l1+l2+l3))).^2+2./P1(2:(l1+l2+l3)).*P1_b(2:(l1+l2+l3))-b./P1(2:(l1+l2+l3)).^2.*P1_b(2:(l1+l2+l3)).^2+b./P1(2:(l1+l2+l3)).*P1_b2(2:(l1+l2+l3)))-P2(2:(l1+l2+l3)).^b.*((log(P2(2:(l1+l2+l3)))+b./P2(2:(l1+l2+l3)).*P2_b(2:(l1+l2+l3))).^2+2./P2(2:(l1+l2+l3)).*P2_b(2:(l1+l2+l3))-b./P2(2:(l1+l2+l3)).^2.*P2_b(2:(l1+l2+l3)).^2+b./P2(2:(l1+l2+l3)).*P2_b2(2:(l1+l2+l3)));
rou_b_a(1)=b*P1(1).^(b-1).*P1_a(1).*(log(P1(1))+b./P1(1).*P1_b(1))+P1(1).^b.*(1./P1(1).*P1_a(1)-b./P1(1).^2.*P1_a(1).*P1_b(1)+b./P1(1).*P1_b_a(1));
rou_b_a(2:(l1+l2+l3))=(b*P1(2:(l1+l2+l3)).^(b-1).*P1_a(2:(l1+l2+l3)).*(log(P1(2:(l1+l2+l3)))+b./P1(2:(l1+l2+l3)).*P1_b(2:(l1+l2+l3)))+P1(2:(l1+l2+l3)).^b.*(1./P1(2:(l1+l2+l3)).*P1_a(2:(l1+l2+l3))-b./P1(2:(l1+l2+l3)).^2.*P1_a(2:(l1+l2+l3)).*P1_b(2:(l1+l2+l3))+b./P1(2:(l1+l2+l3)).*P1_b_a(2:(l1+l2+l3))))-(b*P2(2:(l1+l2+l3)).^(b-1).*P2_a(2:(l1+l2+l3)).*(log(P2(2:(l1+l2+l3)))+b./P2(2:(l1+l2+l3)).*P2_b(2:(l1+l2+l3)))+P2(2:(l1+l2+l3)).^b.*(1./P2(2:(l1+l2+l3)).*P2_a(2:(l1+l2+l3))-b./P2(2:(l1+l2+l3)).^2.*P2_a(2:(l1+l2+l3)).*P2_b(2:(l1+l2+l3))+b./P2(2:(l1+l2+l3)).*P2_b_a(2:(l1+l2+l3))));
rou_a2(1)=b*((b-1)*P1(1).^(b-2).*P1_a(1).^2+P1(1).^(b-1).*P1_a2(1));
rou_a2(2:(l1+l2+l3))=b*((b-1)*P1(2:(l1+l2+l3)).^(b-2).*P1_a(2:(l1+l2+l3)).^2+P1(2:(l1+l2+l3)).^(b-1).*P1_a2(2:(l1+l2+l3)))-b*((b-1)*P2(2:(l1+l2+l3)).^(b-2).*P2_a(2:(l1+l2+l3)).^2+P2(2:(l1+l2+l3)).^(b-1).*P2_a2(2:(l1+l2+l3)));


I(1,1)=lam/u^p*sum(rou,'omitnan');
I(1,2)=0;
I(1,3)=lam/u^(p-1)*sum(rou_b,'omitnan');
I(1,4)=lam/u^(p-1)*sum(rou_a,'omitnan');
I(2,2)=1/(2*lam^2)*(l1+l2+l3);
I(2,3)=1/2*sum(rou_b.*E_d+rou.*E_d_rou.*rou_b,'omitnan');
I(2,4)=1/2*sum(rou_a.*E_d+rou.*E_d_rou.*rou_a,'omitnan');
I(3,3)=1/2*sum((p-1)./rou.^2.*rou_b.^2+(1-p)./rou.*rou_b2+lam*rou_b2.*E_d+2*lam*E_d_rou.*rou_b.^2+lam*rou.*E_d_rou2.*rou_b.^2+lam*rou.*E_d_rou.*rou_b2,'omitnan');
I(3,4)=1/2*sum((p-1)./rou.^2.*rou_a.*rou_b+(1-p)./rou.*rou_b_a+lam*rou_b_a.*E_d+2*lam*E_d_rou.*rou_a.*rou_b+lam*rou.*E_d_rou2.*rou_a.*rou_b+lam*rou.*E_d_rou.*rou_b_a,'omitnan');
I(4,4)=1/2*sum((p-1)./rou.^2.*rou_a.^2+(1-p)./rou.*rou_a2+lam*rou_a2.*E_d+2*lam*E_d_rou.*rou_a.^2+lam*rou.*E_d_rou2.*rou_a.^2+lam*rou.*E_d_rou.*rou_a2,'omitnan');
I(2,1)=I(1,2);
I(3,1)=I(1,3);
I(4,1)=I(1,4);
I(3,2)=I(2,3);
I(4,2)=I(2,4);
I(4,3)=I(3,4);

I=n*I;
Avar=R0_theta*inv(I)*R0_theta';

Std=Avar^0.5;
