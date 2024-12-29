%%绘制CE模型的思想图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%CE模型
clear
subplot(1,2,1)
%退化为凹函数
x=0:0.0001:1;
y1=3*x.^3;
y2=6*x.^3;
y3=8*x.^3;
y4=12*x.^3;


% plot(x,y1)
% hold on
% plot(x,y2)
% plot(x,y3)
% plot(x,y4)

tau1=0.5;
tau2=0.7;
tau3=0.9;
tau4=1;

u1=(3*tau1^3/6)^(1/3);
u2=(6*(u1+(tau2-tau1))^3/8)^(1/3);
u3=(8*(u2+(tau3-tau2))^3/12)^(1/3);

plot(x(1:floor(tau1/0.0001)),y1(1:floor(tau1/0.0001)),'k-','LineWidth',1.5)
hold on
plot(x(floor(tau1/0.0001):floor(tau4/0.0001)),y1(floor(tau1/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u1/0.0001)),y2(1:floor(u1/0.0001)),'k--')
plot(x(floor(u1/0.0001):floor((u1+tau2-tau1)/0.0001)),y2(floor(u1/0.0001):floor((u1+tau2-tau1)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u1+tau2-tau1)/0.0001):floor(tau4/0.0001)),y2(floor((u1+tau2-tau1)/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u2/0.0001)),y3(1:floor(u2/0.0001)),'k--')
plot(x(floor(u2/0.0001):floor((u2+tau3-tau2)/0.0001)),y3(floor(u2/0.0001):floor((u2+tau3-tau2)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u2+tau3-tau2)/0.0001):floor(tau4/0.0001)),y3(floor((u2+tau3-tau2)/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u3/0.0001)),y4(1:floor(u3/0.0001)),'k--')
plot(x(floor(u3/0.0001):floor((u3+tau4-tau3)/0.0001)),y4(floor(u3/0.0001):floor((u3+tau4-tau3)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u3+tau4-tau3)/0.0001):floor(tau4/0.0001)),y4(floor((u3+tau4-tau3)/0.0001):floor(tau4/0.0001)),'k--')


%添加等横线
plot(0:0.01:1,3*0.5^3*ones(length(0:0.01:1)),'k:')
plot(0:0.01:1,6*(u1+0.2)^3*ones(length(0:0.01:1)),'k:')
plot(0:0.01:1,8*(u2+0.2)^3*ones(length(0:0.01:1)),'k:')
%添加等竖线
% plot(u1*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')
% plot(u2*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')
% plot(u3*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')

xlabel('testing time')
ylabel('cumulative degradation')
axis([0 1 0 12])
set(gca,'XTick',0:0.1:1)
set(gca,'YTick',0:1:12)


plot(x(1:floor(tau1/0.0001)),y1(1:floor(tau1/0.0001)),'r-','LineWidth',1.5)
hold on
plot(x(floor(tau1/0.0001):floor(tau2/0.0001)),y2(floor((u1+0.0001)/0.0001):floor((u1+tau2-tau1)/0.0001)),'r-','LineWidth',1.5)
plot(x(floor(tau2/0.0001):floor(tau3/0.0001)),y3(floor((u2-0.0001)/0.0001):floor((u2+tau3-tau2)/0.0001)),'r-','LineWidth',1.5)
plot(x(floor(tau3/0.0001):floor(tau4/0.0001)),y4(floor((u3)/0.0001):floor((u3+tau4-tau3)/0.0001)),'r-','LineWidth',1.5)

annotation('arrow',[0.31,0.34],[0.17,0.17],'Color', 'blue', 'LineWidth', 1);
annotation('arrow',[0.35,0.39],[0.25,0.25],'Color', 'blue', 'LineWidth', 1);
annotation('arrow',[0.37,0.45],[0.4,0.4],'Color', 'blue', 'LineWidth', 1);



subplot(1,2,2)
%退化为凸函数
x=0:0.0001:1;
y1=3*x.^(1/3);
y2=6*x.^(1/3);
y3=8*x.^(1/3);
y4=12*x.^(1/3);


% plot(x,y1)
% hold on
% plot(x,y2)
% plot(x,y3)
% plot(x,y4)

tau1=0.25;
tau2=0.5;
tau3=0.75;
tau4=1.0;

u1=(3*tau1^(1/3)/6)^(3);
u2=(6*(u1+(tau2-tau1))^(1/3)/8)^(3);
u3=(8*(u2+(tau3-tau2))^(1/3)/12)^(3);

plot(x(1:floor(tau1/0.0001)),y1(1:floor(tau1/0.0001)),'k-','LineWidth',1.5)
hold on
plot(x(floor(tau1/0.0001):floor(tau4/0.0001)),y1(floor(tau1/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u1/0.0001)),y2(1:floor(u1/0.0001)),'k--')
plot(x(floor(u1/0.0001):floor((u1+tau2-tau1)/0.0001)),y2(floor(u1/0.0001):floor((u1+tau2-tau1)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u1+tau2-tau1)/0.0001):floor(tau4/0.0001)),y2(floor((u1+tau2-tau1)/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u2/0.0001)),y3(1:floor(u2/0.0001)),'k--')
plot(x(floor(u2/0.0001):floor((u2+tau3-tau2)/0.0001)),y3(floor(u2/0.0001):floor((u2+tau3-tau2)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u2+tau3-tau2)/0.0001):floor(tau4/0.0001)),y3(floor((u2+tau3-tau2)/0.0001):floor(tau4/0.0001)),'k--')


plot(x(1:floor(u3/0.0001)),y4(1:floor(u3/0.0001)),'k--')
plot(x(floor(u3/0.0001):floor((u3+tau4-tau3)/0.0001)),y4(floor(u3/0.0001):floor((u3+tau4-tau3)/0.0001)),'k-','LineWidth',1.5)
plot(x(floor((u3+tau4-tau3)/0.0001):floor(tau4/0.0001)),y4(floor((u3+tau4-tau3)/0.0001):floor(tau4/0.0001)),'k--')


%添加等横线
plot(0:0.01:1,3*tau1^(1/3)*ones(length(0:0.01:1)),'k:')
plot(0:0.01:1,6*(u1+(tau2-tau1))^(1/3)*ones(length(0:0.01:1)),'k:')
plot(0:0.01:1,8*(u2+(tau3-tau2))^(1/3)*ones(length(0:0.01:1)),'k:')
%添加等竖线
% plot(u1*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')
% plot(u2*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')
% plot(u3*ones(length(0:0.01:1)),linspace(0,12,length(0:0.01:1)),'k:')

xlabel('testing time')
ylabel('cumulative degradation')
axis([0 1 0 12])
set(gca,'XTick',0:0.1:1)
set(gca,'YTick',0:1:12)



plot(x(1:floor(tau1/0.0001)),y1(1:floor(tau1/0.0001)),'r-','LineWidth',1.5)
hold on
plot(x(floor(tau1/0.0001):floor(tau2/0.0001)),y2(floor((u1)/0.0001):floor((u1+tau2-tau1)/0.0001)),'r-','LineWidth',1.5)
plot(x(floor(tau2/0.0001):floor(tau3/0.0001)),y3(floor((u2)/0.0001):floor((u2+tau3-tau2)/0.0001)),'r-','LineWidth',1.5)
plot(x(floor(tau3/0.0001):floor(tau4/0.0001)),y4(floor((u3)/0.0001):floor((u3+tau4-tau3)/0.0001)),'r-','LineWidth',1.5)

annotation('arrow',[0.62,0.68],[0.31,0.31],'Color', 'blue', 'LineWidth', 1);
annotation('arrow',[0.65,0.77],[0.44,0.44],'Color', 'blue', 'LineWidth', 1);
annotation('arrow',[0.65,0.85],[0.6,0.6],'Color', 'blue', 'LineWidth', 1);