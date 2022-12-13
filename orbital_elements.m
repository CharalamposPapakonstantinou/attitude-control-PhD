%%
clear all;clc;close all
load('topo.mat','topo','topomap1');
whos topo topomap1




I=[1000 0 0]'; J=[0 1000 0]'; K=[0 0 1000]';


% Eccentricity
Ra=6371+20000; %km
Rp=6371+2000; %km

a=(Ra+Rp)/2;
e=(Ra-Rp)/(Ra+Rp);

v=0:0.01:2*pi;
R=a*(1-e^2)./(1+e*cos(v));
x=R.*cos(v);
y=R.*sin(v);


% Inclination
i=deg2rad(45);

crd=rotx(i)*[x' y' zeros(length(x),1)]';
crd=crd';



% RAAN
RAAN=deg2rad(30);

crd2=rotz(RAAN)*crd';
crd2=crd2';
RAANvector=rotz(RAAN)*I;
% quiver3(0,0,0,10*RAANvector(1),10*RAANvector(2),10*RAANvector(3),'linewidth',5);



% Argument of perigee
w=deg2rad(30);

lq=cross(crd2(1,:),crd2(2,:));
q=quatnormalize([cos(w) sin(w)*lq/norm(lq)]);
crd3=quatmultiply(q,[zeros(length(crd2),1) crd2]);
crd3=crd3(:,2:end);




close all;
plotearth();
plot(x,y) % or polar(v,R)
patch(x,y,'k','facealpha',0.3)
hold on
quiver3(0,0,0,3000,0,0,'linewidth',5);
quiver3(0,0,0,0,3000,0,'linewidth',5);
quiver3(0,0,0,0,0,3000,'linewidth',5);
xlabel('X');ylabel('Y');zlabel('Z');
axis equal
hold on;

patch(crd(:,1),crd(:,2),crd(:,3),'red','facealpha',0.2)
patch(crd2(:,1),crd2(:,2),crd2(:,3),'green','facealpha',0.2,'DiffuseStrength',1)
patch(crd3(:,1)+0.3,crd3(:,2)+0.3,crd3(:,3)+0.3,'blue','facealpha',0.2) %vazw to +0.3 apla gia na mhn exw thema me ta xrwmata sto plotarisma (epeidh ayto to plot einai sto idio plane me to prohgoymeno). kanonika den prepei na yparxei to +0.3




VP=[Rp 0 0]; %perigee vector
VP1=rotx(i)*VP';
quiver3(0,0,0,VP1(1),VP1(2),VP1(3),'r','linewidth',1,'autoscale','off');
plot3(VP1(1),VP1(2),VP1(3),'.red','Markersize',30);
VP2=rotz(RAAN)*VP1;
quiver3(0,0,0,VP2(1),VP2(2),VP2(3),'g','linewidth',1,'autoscale','off');
plot3(VP2(1),VP2(2),VP2(3),'.green','Markersize',30);

q=quatnormalize([cos(w) sin(w)*lq/norm(lq)]);
VP3=quatmultiply(q,[0 VP2']);
VP3=VP3(2:end);
quiver3(0,0,0,VP3(1),VP3(2),VP3(3),'blue','linewidth',1,'autoscale','off');
plot3(VP3(1),VP3(2),VP3(3),'.blue','Markersize',30);





for angl=0:45:359
    hold on
    v=deg2rad(angl);
    Rpos=a*(1-e^2)./(1+e*cos(v));
    POS=[Rpos.*cos(v) Rpos.*sin(v)];
    VP1=rotx(i)*[POS 0]';
    VP2=rotz(RAAN)*VP1;
    q=quatnormalize([cos(w) sin(w)*lq/norm(lq)]);
    VP3=quatmultiply(q,[0 VP2']);
    VP3=VP3(2:end);
    quiver3(0,0,0,VP3(1),VP3(2),VP3(3),'blue','linewidth',1,'autoscale','off');
    plot3(VP3(1),VP3(2),VP3(3),'.blue','Markersize',30);
    if angl==0
        text(VP3(1),VP3(2),VP3(3),'    Perigee','fontweight','bold')
    elseif angl==180
        text(VP3(1),VP3(2),VP3(3),'    Apogee','fontweight','bold')
    else
        text(VP3(1),VP3(2),VP3(3),['    True Anomally = ',num2str(angl),'^o'])
    end
    pause(0.01);
end


%% create gif
filename='orbitelements.gif'
ff=figure(1);
for tt=20:2:380
    view(tt,30)
    pause(0.002);
    frame = getframe(ff);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if tt == 20
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end


%%
figure,
geoshow('landareas.shp','FaceColor',[0.5 1.0 0.5])
xlabel('Longtitude');ylabel('Latitude');

hold on

LL=ecef2lla(crd3*1000);  % *1000 giati thelei metra kai to crd3 einai se km
LL1=LL(:,1);
LL2=LL(:,2)+180; % gia na kanei match me ton xarth poy ksekinaei apo -180

plot(LL2,LL1,'.k','Markersize',4)
plot(LL2-360,LL1,'.k','Markersize',4)
xaxis([-180,180])


hold on

if LL2(1)>=180
    plot(LL2(1)-180,LL1(round(end/2)),'.r','Markersize',30)
    text(LL2(1)-180,LL1(round(end/2)),'   Apogee')
    
    plot(LL2(1)-360,LL1(1),'.b','Markersize',30)
    text(LL2(1)-360,LL1(1),'   Perigee')
elseif LL(1,2)<-180
    plot(LL2(1)+180,LL1(round(end/2)),'.r','Markersize',30)
    text(LL2(1)+180,LL1(round(end/2)),'   Apogee')
    
    plot(LL2(1)+360,LL1(1),'.b','Markersize',30)
    text(LL2(1)+360,LL1(1),'   Perigee')
else
    plot(LL2(1),LL1(round(end/2)),'.r','Markersize',30)
    text(LL2(1),LL1(round(end/2)),'   Apogee')
    
    plot(LL2(1)-180,LL1(1),'.b','Markersize',30)
    text(LL2(1)-180,LL1(1),'   Perigee')
end
% 
% if LL2(round(end/2))>180
%     plot(LL2(round(end/2))-180,LL1(round(end/2)),'.r','Markersize',30)
%     text(LL2(round(end/2))-180,LL1(round(end/2)),'  Perigee')
% elseif LL2(round(end/2))<-180
%     plot(LL2(round(end/2))+180,LL1(round(end/2)),'.r','Markersize',30)
%     text(LL2(round(end/2))+180,LL1(round(end/2)),'  Perigee')
% else
%     plot(LL2(round(end/2)),LL1(round(end/2)),'.r','Markersize',30)
%     text(LL2(round(end/2)),LL1(round(end/2)),'  Perigee')
% end





% for i=1:1
% LLA(i,:) = eci2lla([crd3(1,:)],[2022 1 1 4 52 12.4]);
% end
% plot(LLA(:,2),LLA(:,1),'*r','Markersize',5);













