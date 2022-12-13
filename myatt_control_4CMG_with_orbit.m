clc;clear all;close all;
load('topo.mat','topo','topomap1');
whos topo topomap1
Re = 5;
load('topo.mat');
filename='animated.gif'

[xsphere,ysphere,zsphere] = sphere(50);

w=[0 0 0]';
J=[1 0 0; 0 1 0;0 0 1];
Td=1*[-4 -6 6]';
Tc=[0 0 0]';
w_prev=w;
dt=0.1;u=[0 0 0];
ini=[0 0 -10];
q=[0 ini]'; qref = [0.7071 -0.7071 0 0 ]%[0 ini]; % qref=[0.724 0.427 0.383 0.383];
q_prev=quatnormalize(q');

b=deg2rad(54.73); % rad2deg(acos(1/sqrt(3)))=54.73
bb=pi/2-b;
d1=deg2rad(-90);d2=0;d3=deg2rad(90);d4=0;
dprev=[d1 d2 d3 d4]';
Ts=0.01;
dd=[];TT=[];

vx=[-2 -1 -2 -1 -2 -1 -2 -1];
vy=[-2 -2 -1 -1 -2 -2 -1 -1];
vz=[-2 -2 -2 -2 -1 -1 -1 -1];
H=[vx;vy ;vz]; %Vertices of the cube
S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube
figure(1)
hold on
for i=1:size(S,1)
    Si=S(i,:);
    fill3(H(1,Si),H(2,Si),H(3,Si),[i/size(S,1) 0 1],'facealpha',0.6)
end
axis equal, axis on, hold off, view(20,10)
vnew=[];
stop=20;Integral=[0 0 0]';axiserrorprev=[0;0;0];

%%%%%%%%%%%%%%%%
Kp=20;
Ki=0.0001;
Kwmega=5;

%%%%%%%%%%%%%%%%

endtime=50;
close all

figure('Renderer','painters','Position',[0 900 900 700])
SB1=subplot(2,2,1)
SB2=subplot(2,2,2)
SB3=subplot(2,2,3)
SB4=subplot(2,2,4)

NN=10;
satp=[0:(NN/endtime):NN];
r_cylinder=2;
OPAC=0.9;

for dur=1:endtime
    
    k1=7*cos(satp(dur));
    k2=10*sin(satp(dur))+3;
    k3=0;
    knew=rotx(20)*[k1 k2 k3]';
    vx=[-1 1 -1 1 -1 1 -1 1];
    vy=[-1 -1 1 1 -1 -1 1 1];
    vz=[-1 -1 -1 -1 1 1 1 1];
    H=[vx;vy;vz];
    if dur>stop
        Td=[0 0 0]';
    end
    % expressed in the body frame
    
    wdot=inv(J)*(-skew(w)*J*w+Td+Tc); % the same as wdot=inv(J)*(-cross(w,(J*w)) + Td + Tc);
    w=wdot*dt+w_prev;
    qd=-0.5*quatmultiply([0;w]',q');
    q=qd'*dt+q_prev';
    
    qerror=quatmultiply(quatnormalize(qref),quatnormalize(quatconj(q')))';
    th=2*acos(qerror(1));
    axiserror=qerror(2:end);
    
    if qerror(1)<0
        axiserror=-1.0*axiserror;
    end
    
    Integral=Integral+axiserror;
    if dur==1
        axiserrorprev=axiserror;
    end
    
    u=-Kp*axiserror-Ki*Integral-Kwmega*w;
    
    A=[-cos(b)*cos(d1)   sin(d2)             cos(b)*cos(d3)   -sin(d4);
        -sin(d1)         -cos(b)*cos(d2)     sin(d3)          cos(b)*cos(d4);
        sin(b)*cos(d1)   sin(b)*cos(d2)      sin(b)*cos(d3)   sin(b)*cos(d4)];
    
    h=[-cos(b)*sin(d1)   -cos(d2)             cos(b)*sin(d3)   cos(d4);
        cos(d1)          -cos(b)*sin(d2)      -cos(d3)         cos(b)*sin(d4);
        sin(b)*sin(d1)   sin(b)*sin(d2)       sin(b)*sin(d3)   sin(b)*sin(d4)];
    % NULL SPACE / LCI
    %   sigma=svd(A);
    %   manip(dur)=sigma(3)/sigma(1);
    %   NSgain=?;
    % NULL SPACE / Manipulability
    manip(dur)=sqrt(det(A*A'));
    NSgain=15;
    %%%%%%%
    
    hx=sum(h(1,:));
    hy=sum(h(2,:));
    hz=sum(h(3,:));
    Hvector=[hx hy hz]';
    Hdot=-u-cross(w,Hvector);
    
    ALLCOMBIN=allcomb([d1 d1+0.1 d1-0.1],[d2 d2+0.1 d2-0.1],[d3 d3+0.1 d3-0.1],[d4 d4+0.1 d4-0.1]);
    for i=1:length(ALLCOMBIN)
        d1next=ALLCOMBIN(i,1);d2next=ALLCOMBIN(i,2);d3next=ALLCOMBIN(i,3);d4next=ALLCOMBIN(i,4);
        Anext=[-cos(b)*cos(d1next)   sin(d2next)             cos(b)*cos(d3next)   -sin(d4next);
            -sin(d1next)         -cos(b)*cos(d2next)     sin(d3next)          cos(b)*cos(d4next);
            sin(b)*cos(d1next)   sin(b)*cos(d2next)      sin(b)*cos(d3next)   sin(b)*cos(d4next)];
        manext(i)=sqrt(det(Anext*Anext'));
    end
    [maxw,maxw_index]=max(manext);
    delta_best=ALLCOMBIN(maxw_index,:);
    diafd=delta_best-[d1 d2 d3 d4];
    for i=1:4 % NUMBER OF CMGS
        if diafd(i)~=0
            grad(i)=( maxw-sqrt(det(A*A')) )/diafd(i);
        else
            grad(i)=0;
        end
    end
    d_dot=pinv(A)*Hdot+NSgain*(eye(4,4)-pinv(A)*A)*grad'; % !!! BE CAREFUL FOR -torque !!!
    d=Ts*d_dot+dprev;
    d1=d(1);d2=d(2);d3=d(3);d4=d(4);
    dprev=d;
    
    g1=[sin(b) 0 cos(b)]';
    g2=[0 sin(b) cos(b)]';
    g3=[-sin(b) 0 cos(b)]';
    g4=[0 -sin(b) cos(b)]';
    
    d1_vect=d_dot(1)*g1';
    d2_vect=d_dot(2)*g2';
    d3_vect=d_dot(3)*g3';
    d4_vect=d_dot(4)*g4';
    
    Tc=[cross(h(:,1),d1_vect)+cross(h(:,2),d2_vect)+cross(h(:,3),d3_vect)+cross(h(:,4),d4_vect)]';
    
    vt1=cross(h(:,1),d1_vect);
    vt2=cross(h(:,2),d2_vect);
    vt3=cross(h(:,3),d3_vect);
    vt4=cross(h(:,4),d4_vect);
    
    q_prev=q';
    w_prev=w;
    dd=[dd d_dot];
    TT=[TT Tc];
    
    
    axes(SB1)
    vnew=quatrotate(quatnormalize(q'),ini);
    pt = [0 0 0];
    dir = vnew;
    h = quiver3(pt(1),pt(2),pt(3), dir(1),dir(2),dir(3),'green');
    pt = [0 0 0];
    dir = axiserror;
    h = quiver3(pt(1),pt(2),pt(3), dir(1),dir(2),dir(3),'black');
    quat=q';
    for i=1:length(vx)
        vnew(i,:)=quatrotate(quat,[vx(i) vy(i) vz(i)]); % It uses the dcm matrix
    end
    H=vnew'+knew; %Vertices of the cube
    kk(:,dur)=knew;
    S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube
    
    ff=figure(1)
    s = surf(Re*xsphere,Re*ysphere,Re*zsphere); % create a sphere
    s.CData = topo;  % set color data to topographic data
    s.FaceAlpha=1;
    s.FaceColor = 'texturemap';    % use texture mapping
    s.EdgeColor = 'none';          % remove edges
    s.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces
    s.SpecularStrength = 0.5;      % change the strength of the reflected light
    box on;grid off;
    if dur<stop
        hold on
        for i=1:size(S,1)
            Si=S(i,:);
            fill3(H(1,Si),H(2,Si),H(3,Si),[0 i/size(S,1) 0],'facealpha',1)
        end
        axis([-1 1 -1 1 -1 1]*14)
        axis on,hold off
    else
        hold on
        for i=1:size(S,1)
            Si=S(i,:);
            fill3(H(1,Si),H(2,Si),H(3,Si),[i/size(S,1) 0 0],'facealpha',1)
        end
        axis([-1 1 -1 1 -1 1]*14)
        axis on,hold off,
    end
    view(45-dur/4, 5+dur/7.0);
    xlabel('x');ylabel('y');zlabel('z');
    pause(0.0000005);
    
    
    axes(SB2);
    hold on
    plot(dur,axiserror(3),'.red');
    hold on
    plot(dur,axiserror(2),'.green');
    hold on
    plot(dur,axiserror(1),'.blue');
    if dur==endtime legend('error in axis z', 'error in axis y', 'error in axis x','Location','southeast') ,  end
    axis([0 endtime -1 1]);
    title('Error in x,y,z')
    xlabel('time')
    
    
    axes(SB3);clc
    hold on
    plot(dur,rad2deg(d1),'.red');
    hold on
    plot(dur,rad2deg(d2),'.green');
    hold on
    plot(dur,rad2deg(d3),'.blue');
    hold on
    plot(dur,rad2deg(d4),'.black');
    if dur==endtime legend('d1', 'd2', 'd3','d4','Location','northeast') ,end
    axis([0 endtime -250 200]);
    title('Gimbal angles ')
    xlabel('time')
    ylabel('theta(degrees)')
    
    %%%  Visualization %%%%
    
    axes(SB4)
    hold(SB4,'off');
   %%
   
%    f1= figure('Renderer','painters','Units','inches','Position',[0 0 3.5 2.5], 'PaperPositionMode','auto','Name','MYFIG1');

    view(45,80)
    
    %%%%% CMG 1 %%%%%
    [X,Y,Z] = cylinder(r_cylinder);
    
    CMG3x=reshape(X,[1 2*length(X)]);
    CMG3y=reshape(Y,[1 2*length(Y)]);
    CMG3z=reshape(Z-0.5,[1 2*length(Z)]);
    CMG3=rotz(pi)*roty(bb)*rotx(pi/2-d1)*[CMG3x;CMG3y;CMG3z]+[4 0 0]';
    
    CMG3xnew=reshape(CMG3(1,:),[2 length(CMG3(1,:))/2]);
    CMG3ynew=reshape(CMG3(2,:),[2 length(CMG3(2,:))/2]);
    CMG3znew=reshape(CMG3(3,:),[2 length(CMG3(3,:))/2]);
    
    surf(CMG3xnew,CMG3ynew,CMG3znew,'facecolor','r','LineStyle','none','facealpha',OPAC);axis auto;
    hold on
    fill3(CMG3xnew(1,:),CMG3ynew(1,:),CMG3znew(1,:),'r','facealpha',OPAC)
    fill3(CMG3xnew(2,:),CMG3ynew(2,:),CMG3znew(2,:),'r','facealpha',OPAC)
    
    xl1=  [0 0;0 0;-1 -1;-1 1;-1 1;-1 -1;1 1;-1 1];
    yl1=1.5*[-3 -2;2 3;0 0;-2 -2;2 2;-2 2;-2 2;0 0];
    zl1=zeros(length(xl1),2);

    
    xyzlineA=rotz(pi)*roty(bb)*rotx(pi/2-d1)*rotx(pi/2)*rotz(pi/2)*[xl1(:,1) yl1(:,1) zl1(:,1)]'+[4 0 0]';
    xyzlineB=rotz(pi)*roty(bb)*rotx(pi/2-d1)*rotx(pi/2)*rotz(pi/2)*[xl1(:,2) yl1(:,2) zl1(:,2)]'+[4 0 0]';
    xl1new=[xyzlineA(1,:);xyzlineB(1,:)];
    yl1new=[xyzlineA(2,:);xyzlineB(2,:)];
    zl1new=[xyzlineA(3,:);xyzlineB(3,:)];
    for i=1:length(xl1new)
        line(xl1new(:,i),yl1new(:,i),zl1new(:,i),'linewidth',4);
        hold on
    end
    text(6,0,4,strcat('\delta_1=',num2str(rad2deg(d1),3),char(176)),'fontsize',15)
    
%     quiver3(4,0,0,vt1(1),vt1(2),vt1(3),'LineWidth',2,'Color',[0 1 0],'AutoScale','on','AutoScaleFactor',1.5)
% 
%     cross(h(:,1),d1_vect)
%     vt1=vt1(1)-4;
%     quiver3(0,0,0,g1(1),g1(2),g1(3),'LineWidth',3,'Color',[1 0 0]);hold on
%     line([0 vt3(1)],[0 vt3(2)],[0 vt3(3)])
    
    %%%%% CMG 2 %%%%%
    [X,Y,Z] = cylinder(r_cylinder);
    
    CMG4x=reshape(X,[1 2*length(X)]);
    CMG4y=reshape(Y,[1 2*length(Y)]);
    CMG4z=reshape(Z-0.5,[1 2*length(Z)]);
    CMG4=rotz(3*pi/2)*roty(bb)*rotx(pi/2-d2)*[CMG4x;CMG4y;CMG4z]+[0 4 0]';
    
    CMG4xnew=reshape(CMG4(1,:),[2 length(CMG4(1,:))/2]);
    CMG4ynew=reshape(CMG4(2,:),[2 length(CMG4(2,:))/2]);
    CMG4znew=reshape(CMG4(3,:),[2 length(CMG4(3,:))/2]);
    
    surf(CMG4xnew,CMG4ynew,CMG4znew,'facecolor','r','LineStyle','none','facealpha',OPAC);axis auto;
    hold on
    fill3(CMG4xnew(1,:),CMG4ynew(1,:),CMG4znew(1,:),'r','facealpha',OPAC)
    fill3(CMG4xnew(2,:),CMG4ynew(2,:),CMG4znew(2,:),'r','facealpha',OPAC)
    
    xl1=  [0 0;0 0;-1 -1;-1 1;-1 1;-1 -1;1 1;-1 1];
    yl1=1.5*[-3 -2;2 3;0 0;-2 -2;2 2;-2 2;-2 2;0 0];
    zl1=zeros(length(xl1),2);
    
    
    xyzlineA=roty(pi/2)*rotx(pi/2)*roty(bb)*rotz(pi/2-d2)*rotx(pi/2)*[xl1(:,1) yl1(:,1) zl1(:,1)]'+[0 4 0]';
    xyzlineB=roty(pi/2)*rotx(pi/2)*roty(bb)*rotz(pi/2-d2)*rotx(pi/2)*[xl1(:,2) yl1(:,2) zl1(:,2)]'+[0 4 0]';
    xl1new=[xyzlineA(1,:);xyzlineB(1,:)];
    yl1new=[xyzlineA(2,:);xyzlineB(2,:)];
    zl1new=[xyzlineA(3,:);xyzlineB(3,:)];
    for i=1:length(xl1new)
        line(xl1new(:,i),yl1new(:,i),zl1new(:,i),'linewidth',4);
        hold on
    end
 
    text(0,6,4,strcat('\delta_2=',num2str(rad2deg(d2),3),char(176)),'fontsize',15)
%     quiver3(0,4,0,vt2(1),vt2(2),vt2(3),'LineWidth',2,'Color',[0 1 0],'AutoScale','on','AutoScaleFactor',1.5)

    
    %%%%% CMG3 %%%%%%%%%%%%%
    [X,Y,Z] = cylinder(r_cylinder);
    %     h=mesh(X,Y,Z,'facecolor',[1 0 0])
    
    CMG1x=reshape(X,[1 2*length(X)]);
    CMG1y=reshape(Y,[1 2*length(Y)]);
    CMG1z=reshape(Z-0.5,[1 2*length(Z)]);
    CMG1=roty(bb)*rotx(pi/2-d3)*[CMG1x;CMG1y;CMG1z]+[-4 0 0]';
    
    CMG1xnew=reshape(CMG1(1,:),[2 length(CMG1(1,:))/2]);
    CMG1ynew=reshape(CMG1(2,:),[2 length(CMG1(2,:))/2]);
    CMG1znew=reshape(CMG1(3,:),[2 length(CMG1(3,:))/2]);
    
    surf(CMG1xnew,CMG1ynew,CMG1znew,'facecolor','r','LineStyle','none','facealpha',OPAC);axis auto;
    hold on
    fill3(CMG1xnew(1,:),CMG1ynew(1,:),CMG1znew(1,:),'r','facealpha',OPAC)
    fill3(CMG1xnew(2,:),CMG1ynew(2,:),CMG1znew(2,:),'r','facealpha',OPAC)
    
    xl1=  [0 0;0 0;-1 -1;-1 1;-1 1;-1 -1;1 1;-1 1];
    yl1=1.5*[-3 -2;2 3;0 0;-2 -2;2 2;-2 2;-2 2;0 0];
    zl1=zeros(length(xl1),2);
    
    xyzlineA=roty(bb)*rotx(pi/2-d3)*rotx(pi/2)*rotz(pi/2)*[xl1(:,1) yl1(:,1) zl1(:,1)]'+[-4 0 0]';
    xyzlineB=roty(bb)*rotx(pi/2-d3)*rotx(pi/2)*rotz(pi/2)*[xl1(:,2) yl1(:,2) zl1(:,2)]'+[-4 0 0]';
    xl1new=[xyzlineA(1,:);xyzlineB(1,:)];
    yl1new=[xyzlineA(2,:);xyzlineB(2,:)];
    zl1new=[xyzlineA(3,:);xyzlineB(3,:)];
    for i=1:length(xl1new)
        line(xl1new(:,i),yl1new(:,i),zl1new(:,i),'linewidth',4);
        hold on
    end
    axis equal
    hold on
    text(-6,0,-3,strcat('\delta_3=',num2str(rad2deg(d3),3),char(176)),'fontsize',15)
%     quiver3(-4,0,0,vt3(1),vt3(2),vt3(3),'LineWidth',2,'Color',[0 1 0],'AutoScale','on','AutoScaleFactor',1.5)

%     cross(h(:,1),d1_vect)
%     vt1=vt1(1)-4;
%     quiver3(0,0,0,g3(1),g3(2),g3(3),'LineWidth',3,'Color',[1 0 0]);hold on
%     line([0 vt3(1)],[0 vt3(2)],[0 vt3(3)])
    
    
    %%%%% CMG 4 %%%%%
    [X,Y,Z] = cylinder(r_cylinder);
    
    CMG2x=reshape(X,[1 2*length(X)]);
    CMG2y=reshape(Y,[1 2*length(Y)]);
    CMG2z=reshape(Z-0.5,[1 2*length(Z)]);
    CMG2=rotz(pi/2)*roty(bb)*rotx(pi/2-d4)*[CMG2x;CMG2y;CMG2z]+[0 -4 0]';
    
    CMG2xnew=reshape(CMG2(1,:),[2 length(CMG2(1,:))/2]);
    CMG2ynew=reshape(CMG2(2,:),[2 length(CMG2(2,:))/2]);
    CMG2znew=reshape(CMG2(3,:),[2 length(CMG2(3,:))/2]);
    
    surf(CMG2xnew,CMG2ynew,CMG2znew,'facecolor','r','LineStyle','none','facealpha',OPAC);axis auto;
    hold on
    fill3(CMG2xnew(1,:),CMG2ynew(1,:),CMG2znew(1,:),'r','facealpha',OPAC)
    fill3(CMG2xnew(2,:),CMG2ynew(2,:),CMG2znew(2,:),'r','facealpha',OPAC)
    
    xl1=  [0 0;0 0;-1 -1;-1 1;-1 1;-1 -1;1 1;-1 1];
    yl1=1.5*[-3 -2;2 3;0 0;-2 -2;2 2;-2 2;-2 2;0 0];
    zl1=zeros(length(xl1),2);
    
    
    xyzlineA=rotx(-bb)*roty(pi/2-d4)*roty(pi/2)*[xl1(:,1) yl1(:,1) zl1(:,1)]'+[0 -4 0]';
    xyzlineB=rotx(-bb)*roty(pi/2-d4)*roty(pi/2)*[xl1(:,2) yl1(:,2) zl1(:,2)]'+[0 -4 0]';
    xl1new=[xyzlineA(1,:);xyzlineB(1,:)];
    yl1new=[xyzlineA(2,:);xyzlineB(2,:)];
    zl1new=[xyzlineA(3,:);xyzlineB(3,:)];
    for i=1:length(xl1new)
        line(xl1new(:,i),yl1new(:,i),zl1new(:,i),'linewidth',4);
        hold on
    end
    text(0,-6,-3,strcat('\delta_4=',num2str(rad2deg(d4),3),char(176)),'fontsize',15)
    
%     quiver3(0,-4,0,vt4(1),vt4(2),vt4(3),'LineWidth',2,'Color',[0 1 0],'AutoScale','on','AutoScaleFactor',1.5)

    
    
    title('CMGs | Visualization');xlabel('x');ylabel('y');zlabel('z');
    axis ([-8 8 -8 8 -5 5]);
    %%
    
    

    frame = getframe(ff);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if dur == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

    
end
axes(SB1)
hold on
plot3(SB1,kk(1,:),kk(2,:),kk(3,:),'.red');

figure('Renderer','painters','Position',[00 0 900 350])
subplot(1,3,1)
plot([0:endtime-1],dd,'-')
ylabel('CMG Gimbal Rate')
xlabel('time')
legend('CMG_1', 'CMG_2', 'CMG_3','CMG_4','Location','northeast')

subplot(1,3,2)
plot([0:endtime-1],TT,'-')
ylabel('CMG Torque')
xlabel('time')
legend('x', 'y', 'z','Location','northeast')

subplot(1,3,3)
plot([0:endtime-1],manip,'-')
ylabel('Manipulability')
xlabel('time')




%%
