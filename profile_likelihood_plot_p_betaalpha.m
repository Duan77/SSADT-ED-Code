clear
%似然函数关于p的截面图
p1=-2:0.01:-0.01;
p2=1.01:0.01:1.99;%%0 1 2三点处比较特殊，是特例
p3=2.01:0.01:3;
for i=1:length(p1)
    theta=[p1(i),0.484549499028519,1.70641297762755];
    fun1(i)=-likefun_step_stress(theta);
end

hold on
for j=1:length(p2)
    theta=[p2(j),0.484549499028519,1.70641297762755];
    fun2(j)=-likefun_step_stress(theta);
end

for k=1:length(p3)
    theta=[p3(k),0.484549499028519,1.70641297762755];
    fun3(k)=-likefun_step_stress(theta);
end

[maxVal, linearInd] = max(fun2); % 最大值和线性索引
%[row, col] = ind2sub(size(p2), linearInd); % 行和列索引



%似然函数关于beta alpha的截面图

beta=0.01:0.01:2;
alpha=0.01:0.01:4;
[x,y] = meshgrid(beta,alpha);
for i=1:length(beta)
    for j=1:length(alpha)
        theta1=[1.54458138438586,beta(i),alpha(j)];
        z(i,j) = -likefun_step_stress(theta1);
    end
end


for i=1:length(beta)
 for j=1:length(alpha)
   if z(i,j)==max(max(z))
      row=i;
      col=j;
   end
 end
end


subplot(1,2,1)
plot(p1,fun1,'r*','markersize',5)
hold on
plot(p2,fun2,'r*','markersize',5)
plot(p3,fun3,'r*','markersize',5)
plot(p2(linearInd),fun2(linearInd),'k.','markersize',30)
%plot(-2:0.01:p2(linearInd),repelem(maxVal,length(-2:0.01:p2(linearInd))),'k:','LineWidth',2)
%plot(repelem(p2(linearInd),length(-1100:1:fun2(linearInd))),-1100:1:fun2(linearInd),'k:','LineWidth',2)
text(p2(linearInd),fun2(linearInd)+30,'p=1.54','FontSize',10)
%text(p2(linearInd),fun2(linearInd)+15,'lnL=-453.6863','FontSize',10)

xlabel('Parameter p')
ylabel('Proflie log-likelihood')
% axis([-2 3 -1100 -400])
% set(gca,'XTick',-2:0.5:3)
% set(gca,'YTick',-1100:100:-400)



subplot(1,2,2)
% 绘制三维曲面
meshc(x, y, z');

hold on
plot3(beta(row),alpha(col),z(row,col),'.','color','k','markersize',30)
text(beta(row),alpha(col),z(row,col)+280,'\beta=0.48')
text(beta(row),alpha(col),z(row,col)+160,'\alpha=1.71')

%axis square; %使三个坐标轴等长
%grid on%添加网格 off丢掉网格线
colormap jet
colorbar
%shading interp 
xlabel('Parameter \beta')
ylabel('Parameter \alpha')
zlabel('Proflie log-likelihood')
