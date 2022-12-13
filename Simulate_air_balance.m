clc;clear all;close all;

warning off
MYTABLE=[];

w=[0 0 0]';
J=[0.00809 0 0; 0 0.00817365 0.0;0 0 0.01335037]; % Mass=24kg , Side of Cube=0.5m
Td=[0 0 0]';
Tc=[0 0 0]';
w_prev=w;
u=[0 0 0];
ini=[1 0 0];
% q=[1 0 0 0]';
qref=[1 0 0 0];q=qref';
% q=0*eul2quat([0.0041  0.0601 -0.3405 ])';
q_prev=quatnormalize(q');

b=deg2rad(54.73);
d1=deg2rad(-70);
d2=0;
d3=deg2rad(75);
d4=0;
% d1=pi/1.91;d2=d1;d3=d1;d4=d1; % external elliptic
dprev=[d1 d2 d3 d4]';

dd=[];TT=[];d1all=[];d2all=[];d3all=[];d4all=[];axer=[];

vx=[-1 1 -1 1 -1 1 -1 1];
vy=[-1 -1 1 1 -1 -1 1 1];
vz=[-1 -1 -1 -1 1 1 1 1];
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
Integral=[0 0 0]';axiserrorprev=[0;0;0];

%%%%%%%%%%%%%%%%
gfactor=1;
Kp=8/gfactor;
Ki=0*0.00001/gfactor;
Kwmega=15/gfactor;
%%%%%%%%%%%%%%%%
dt=0.01;
endtime=20/dt;
timeinsec=[0:dt:endtime*dt-dt];
stop=0/dt;
plot_motion=0;

close all
%
SB1=figurem('Renderer','painters','Position',[0 100 700 700])
% SB2=subplot(1,3,2)
% SB3=subplot(1,3,3)
clc;
h0=1;
epsilon=deg2rad(1.5);
singpoint=0;
Favoid=[0 0 0]';
dqtot=[0 0 0 0];
tsoption=0;
Favoidprev=[0;0;0];
Hvector=[0;0;0];sc=0;tc=0;sss=0;

cross_sum=0;gravskewtot=0;rightterm2tot=0;LMAT=[];RMAT=[];
PHIMAT=[];DWMAT=[];Tdd=[];
fi=0;theta=0;psi=0;
fiprev=0;thetaprev=0;psiprev=0;
w0prev = w(1); w1prev = w(2); w2prev = w(3);
OPTION=0; % 0 = only r  ,   1 =r and J
tic

for dur=1:endtime
%         att_angle(dur,:)=[psi theta fi];%quat2eul((q'));
%         fi=w(1)*dt+fiprev;
%         theta=w(2)*dt+thetaprev;
%         psi=w(3)*dt+psiprev;
%         
%         if abs(psi)>2*pi, psi=sign(psi)*(mod(psi,2*pi));, end;
%         if abs(theta)>2*pi, theta=sign(theta)*(mod(theta,2*pi));, end;
%         if abs(fi)>2*pi, fi=sign(fi)*(mod(fi,2*pi));, end;        
    
    att_angle(dur,:)=quat2eul((q'));
    fi=att_angle(dur,3);
    theta=att_angle(dur,2);
    psi=att_angle(dur,1);
    
    rx= -0.00023163;ry=0.0002117;rz= -1.7625*1e-3;
    Mass=1.166;
    if dur>stop
        Td=Mass*9.81*[rz*cos(theta)*sin(fi)-ry*cos(theta)*cos(fi);rz*sin(theta)+rx*cos(theta)*cos(fi);-ry*sin(theta)-rx*cos(theta)*sin(fi)];
        Tdd=[Tdd,Td];
        qref = [0.70711 -0.70711 0 0]; %[0.70711 -0.70711 0 0.5],[0.70711 -0.70711 0.5 0.5],[0.70711 -0.70711 -1.5 -1.5]
        %       qref=[0.70711 -0.70711 -1.5 -1.5];
    end
    % expressed in the body frame
    wdot=inv(J)*(-skew(w)*J*w + Td -OPTION*cross(w,Hvector)); % the same as wdot=inv(J)*(-cross(w,(J*w)) + Td + Tc);
    w=wdot*dt+w_prev;
    
    qd=0.5*quatmultiply(q',[0;w]');
    if ( norm(w)~=0 )
        q=quatmultiply(q_prev,[cos(norm(w)*dt/2);(w/norm(w))*sin(norm(w)*dt/2)]')';
    end
    
    qerror=quatmultiply(quatconj(quatnormalize(qref)),quatnormalize((q')))';
    eulererror = quat2eul(qerror');
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
    
    d1=pi/2;
    d2=pi/2;
    d3=pi/2;
    d4=pi/2;
    
    A=[-cos(b)*cos(d1)   sin(d2)             cos(b)*cos(d3)   -sin(d4);
        -sin(d1)         -cos(b)*cos(d2)     sin(d3)          cos(b)*cos(d4);
        sin(b)*cos(d1)   sin(b)*cos(d2)      sin(b)*cos(d3)   sin(b)*cos(d4)];
    
    h=0.001041332*[-cos(b)*sin(d1)   -cos(d2)             cos(b)*sin(d3)   cos(d4);
                    cos(d1)          -cos(b)*sin(d2)      -cos(d3)         cos(b)*sin(d4);
                    sin(b)*sin(d1)   sin(b)*sin(d2)       sin(b)*sin(d3)   sin(b)*sin(d4)];
    
    
    hx=sum(h(1,:));
    hy=sum(h(2,:));
    hz=sum(h(3,:));
    Hvector=[hx hy hz]';
    
    
    
   
    
    if OPTION==1
        
        WMEGA=[w(1), 0, 0, w(2), w(3), 0;
            0, w(2), 0, w(1), 0, w(3);
            0, 0, w(3), 0, w(1), w(2)];
        cross_sum = cross_sum + (skew(w) * WMEGA) * dt;
        leftterm1 = WMEGA + cross_sum;
        gravskew = skew(9.81*[sin(theta);-sin(fi)*cos(theta);-cos(fi)*cos(theta)]);
        gravskewtot = gravskewtot + gravskew * dt;
        leftterm = [leftterm1 , gravskewtot];
        
        rightterm2now = cross(w, Hvector) ;
        rightterm2tot = rightterm2tot + rightterm2now * dt;
        rightterm = -0*Hvector - rightterm2tot;   % set Hvector =0 to work but why??
        ss0(:,dur)=-0*Hvector - rightterm2tot;
        ss1(:,dur)=-Hvector - rightterm2tot;
        
        LMAT=[LMAT;leftterm];
        RMAT=[RMAT;rightterm];
        
        RESMAT(:,dur)=pinv(LMAT)*RMAT;
        
%         fi=w(1)*dt+fiprev;  anti gia to quat2eul
%         theta=w(2)*dt+thetaprev;
%         psi=w(3)*dt+psiprev;

    
        
    elseif OPTION==0
        
        Ixx=J(1,1);Iyy=J(2,2);Izz=J(3,3);
        M=Mass;
        g=9.81;
        
        ph12 = - M * g * dt * ( cos(fi) * cos(theta) + cos(fiprev) * cos(thetaprev) ) / (2 * Ixx);
        
        ph13 = M * g * dt * ( sin(fi) * cos(theta) + sin(fiprev) * cos(thetaprev) ) / (2 * Ixx);
        
        ph21 = M * g * dt * ( cos(fi) * cos(theta) + cos(fiprev) * cos(thetaprev) ) / (2 * Iyy);
        
        ph23 = M * g * dt * ( sin(theta) + sin(thetaprev) ) / (2 * Iyy);
        
        ph31 = - M * g * dt * ( sin(fi) * cos(theta) + sin(fiprev) * cos(thetaprev) ) / (2 * Izz);
        
        ph32 = - M * g * dt * ( sin(theta) + sin(thetaprev) ) / (2 * Izz);
        
        dwx = w(1) - w0prev;
        dwy = w(2) - w1prev;
        dwz = w(3) - w2prev;
        
        w0prev = w(1);
        w1prev = w(2);
        w2prev = w(3);
        
        
        phi = [0, ph12, ph13; ph21, 0, ph23; ph31, ph32, 0];
        dw = [dwx, dwy, dwz]';
        
        PHIMAT=[PHIMAT;phi];
        DWMAT=[DWMAT;dw];
        RES2MAT(:,dur)=pinv(PHIMAT)*DWMAT;
       
%         fi=w(1)*dt+fiprev; anti gia to quat2eul
%         theta=w(2)*dt+thetaprev;
%         psi=w(3)*dt+psiprev;
        
    end
    
    
    
    
    
    fiprev = fi;
    thetaprev = theta;
    psiprev = psi;
    
    
    
    
    
    
    Hdot=[0,0,0]';
    
    manip(dur)=sqrt(det(A*A'));
    
    Ahash=pinv(A);
    d_dot=[0,0,0,0]';
    
    %%%%%%%%%%%% SATURATION %%%%%%%%%%%%%%
    ddmax=max(abs(d_dot));
    ddlim=deg2rad(50);
    ddgamma=ddmax/ddlim;
    if ddgamma>1  d_dot=d_dot/ddgamma; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % !!! BE CAREFUL FOR -torque !!!
    d=[d1 d2 d3 d4]';
    d1=d(1);d2=d(2);d3=d(3);d4=d(4);
    dprev=d;
    
    
    q_prev=q';
    w_prev=w;
    dd=[dd d_dot];
    
    
    
    
    
    
    
    
    
    
    
    
    if plot_motion==1 && mod(dur,24)==0
        axes(SB1)
        vnew=quatrotate(quatconj(q'),ini);
        pt = [0 0 0];
        dir = vnew;
        hq = quiver3(pt(1),pt(2),pt(3), dir(1),dir(2),dir(3),'green');
        pt = [0 0 0];
        dir = axiserror;
        hq = quiver3(pt(1),pt(2),pt(3), dir(1),dir(2),dir(3),'black');
        quat=q';
        for i=1:length(vx)
            vnew(i,:)=quatrotate(quatconj(quat),[vx(i) vy(i) vz(i)]); % It uses the dcm matrix
        end
        H=vnew'; %Vertices of the cube
        S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube
        
        figure(1)
        if dur<stop
            hold on
            for i=1:size(S,1)
                Si=S(i,:);
                fill3(H(1,Si),H(2,Si),H(3,Si),[0 i/size(S,1) 0],'facealpha',1)
            end
            axis([-3 3 -3 3 -3 3]), axis on, hold off, view(20,10)
        else
            hold on
            for i=1:size(S,1)
                Si=S(i,:);
                fill3(H(1,Si),H(2,Si),H(3,Si),[i/size(S,1) 0 0],'facealpha',1)
            end
            axis([-3 3 -3 3 -3 3]), axis on, hold off
        end
        xlabel('x');ylabel('y');zlabel('z'); view(120,45);
        pause(0.000000000000000001);
    end
    
    
    
    d1all=[d1all d1];d2all=[d2all d2];d3all=[d3all d3];d4all=[d4all d4];
    axer=[axer ;axiserror(1) axiserror(2) axiserror(3)];
    eulerror(dur,:)=eulererror';
    %     manip(dur)=sqrt(det(A*A'));
    sigma=svd(A);
    lci_index(dur)=sigma(3)/sigma(1);
    hforplot(:,dur)=Hvector;
    wforplot(:,dur)=w;
    Hdotplot(:,dur)=Hdot;
    
    hh(dur)=norm(Hvector);
    
end
%
if OPTION==1
    plot(RESMAT','Linewidth',2);
    legend('J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}','mrx','mry','mrz')
    mrx=RESMAT(7,end);
    mry=RESMAT(8,end);
    mrz=RESMAT(9,end);
    RESMAT(:,end)
    disp('- - -')
    rx=(mrx/Mass) %in m
    ry=(mry/Mass) %in m
    rz=(mrz/Mass) %in m
    ylim([-0.1,0.1])
elseif OPTION==0
    figure,
    plot(RES2MAT','Linewidth',2);
    ylim([-0.1,0.1])
    legend('mrx','mry','mrz')
    RES2MAT(:,end)
    disp('- - -')
end

tic

% f1=figurem('Renderer','painters','Position',[0 0 1000 1600]);
f1= figure('Renderer','painters','Units','inches','Position',[0 0 7 8], 'PaperPositionMode','auto');
xstyle='-black';ystyle='--black';zstyle=':black';text_X=-0.05;text_Y=-0.14;textsize=9; Nstep=300;

gca=subplot(4,2,1);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,rad2deg(att_angle(:,3)),xstyle);
% plot(timeinsec,rad2deg(att_angle(:,2)),ystyle);
% plot(timeinsec,rad2deg(att_angle(:,1)),zstyle);
xlabel('Time, s')
ylabel('Attitude angle, deg')
legend('x','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=(yLimits(2)-yLimits(1))/Nstep;
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'a)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yticks([-100:20:20])


gca=subplot(4,2,2);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
% plot(timeinsec,rad2deg(axer(:,1)),xstyle);
% plot(timeinsec,rad2deg(axer(:,2)),ystyle);
% plot(timeinsec,rad2deg(axer(:,3)),zstyle);
% plot(timeinsec,rad2deg(att_angle(:,3)),xstyle);
plot(timeinsec,rad2deg(att_angle(:,2)),ystyle);
plot(timeinsec,rad2deg(att_angle(:,1)),zstyle);
ylabel('Attitude Angle, deg')
xlabel('Time, s')
legend('y','z','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'b)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)

gca=subplot(4,2,3);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,rad2deg(wforplot(1,:)),xstyle);
plot(timeinsec,rad2deg(wforplot(2,:)),ystyle);
plot(timeinsec,rad2deg(wforplot(3,:)),zstyle);
xlabel('Time, s')
ylabel('Angular Velocity, deg/s')
legend('x', 'y', 'z','Location','southeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'c)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)


gca=subplot(4,2,4);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,hforplot(1,:),xstyle);
plot(timeinsec,hforplot(2,:),ystyle);
plot(timeinsec,hforplot(3,:),zstyle);
xlabel('Time, s')
ylabel('CMG Momentum, Nms')
legend('x', 'y', 'z','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'d)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)


gca=subplot(4,2,5);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,rad2deg(d1all),xstyle);
plot(timeinsec,rad2deg(d2all),ystyle);
plot(timeinsec,rad2deg(d3all),zstyle);
plot(timeinsec,rad2deg(d4all),'-.black');
ylabel('Gimbal Angles, deg')
xlabel('Time, s')
legend('CMG_1', 'CMG_2', 'CMG_3','CMG_4','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'e)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)


% figure('Renderer','painters','Position',[20 50 1800 300])
gca=subplot(4,2,6);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,rad2deg(dd(1,:)),xstyle)
plot(timeinsec,rad2deg(dd(2,:)),ystyle)
plot(timeinsec,rad2deg(dd(3,:)),zstyle)
plot(timeinsec,rad2deg(dd(4,:)),'-.black','linewidth',0.1)
ylabel('Gimbal Rates, deg/s')
xlabel('Time, s')
legend('CMG_1', 'CMG_2', 'CMG_3','CMG_4','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'f)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yaxis([-50 50])


gca=subplot(4,2,7);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
% plot(timeinsec,TT(1,:),xstyle)
% plot(timeinsec,TT(2,:),ystyle)
% plot(timeinsec,TT(3,:),zstyle)
ylabel('CMG Torque, Nm')
xlabel('Time, s')
legend('x', 'y', 'z','Location','northeast','orientation','vertical')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'g)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yaxis([-1 1])

gca=subplot(4,2,8);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',textsize);
hold on
plot(timeinsec,manip,'-black')
ylabel('Manipulability Index')
xlabel('Time, s')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'h)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)

%%