clc;clear all;close all;
load('surface.mat');
load('satellite_coordinates.mat');
[vertices,faces,normals,name] = stlRead('ascii_sat.stl');%satellite_keychain

w=[0 0 0]';
J=1*[1 0 0; 0 1 0;0 0 1]; % Mass=24kg , Side of Cube=0.5m
Td=0*[0 10 0]';
Tc=[0 0 0]';
w_prev=w;
dt=0.01;u=[0 0 0];
ini=[1 0 0];
% q=[1 0 0 0]'; 
qref=[1 0 0 0];q=qref';

q_prev=quatnormalize(q');

b=deg2rad(54.73);
d1=deg2rad(-70);
d2=0;
d3=deg2rad(75);
d4=0;
dprev=[d1 d2 d3 d4]';

dd=[];TT=[];d1all=[];d2all=[];d3all=[];d4all=[];axer=[];

vx=[-1 1 -1 1 -1 1 -1 1];
vy=[-1 -1 1 1 -1 -1 1 1];
vz=[-1 -1 -1 -1 1 1 1 1];
cstep=1;
vx=vertices(1:cstep:end,1)-28.86;     %% prepei na paw to antikeimeno sto kentro
vy=vertices(1:cstep:end,2)-0.1;
vz=vertices(1:cstep:end,3)-0.11;
H=[vx;vy ;vz]; %Vertices of the cube
figure(1)
hold on
object.vertices = [vx vy vz];
object.faces = [faces(1:cstep:end,1) faces(1:cstep:end,2) faces(1:cstep:end,3)];
patch(object,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
camlight('headlight');
material('dull');
axis(1*[-1 1 -1 1 -1 1]), hold off, view(20,10)

vnew=[];
Integral=[0 0 0]';axiserrorprev=[0;0;0];

%%%%%%%%%%%%%%%%
gfactor=1;
Kp=20/gfactor;
Ki=0.00001/gfactor;
Kwmega=5/gfactor;
%%%%%%%%%%%%%%%%
endtime=20/dt;
timeinsec=[0:dt:endtime*dt-dt];
stop=5/dt;
plot_motion=1;

close all
% 
SB1=figurem('Renderer','painters','Position',[0 100 700 600])
clc;
h0=1;
epsilon=deg2rad(1.5);
singpoint=0;
for dur=1:endtime

    if dur>stop
        Td=[0 0 0]';
        qref = [0.70711 -0.30711 0 0]; % 
    end
    % expressed in the body frame
    wdot=inv(J)*(-skew(w)*J*w + Td + Tc); % the same as wdot=inv(J)*(-cross(w,(J*w)) + Td + Tc);
    w=wdot*dt+w_prev;
    
    qd=0.5*quatmultiply(q',[0;w]');
    if ( w(1)~=0 && w(2)~=0 && w(3)~=0 )
         q=quatmultiply(q_prev,[cos(norm(w)*dt/2);(w/norm(w))*sin(norm(w)*dt/2)]')';
    end
    qerror=quatmultiply(quatconj(quatnormalize(qref)),quatnormalize((q')))';
    
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
    
    h=h0*[-cos(b)*sin(d1)   -cos(d2)             cos(b)*sin(d3)   cos(d4);
        cos(d1)          -cos(b)*sin(d2)      -cos(d3)         cos(b)*sin(d4);
        sin(b)*sin(d1)   sin(b)*sin(d2)       sin(b)*sin(d3)   sin(b)*sin(d4)];
    
    
    hx=sum(h(1,:));
    hy=sum(h(2,:));
    hz=sum(h(3,:));
    Hvector=[hx hy hz]';
    Hdot=-u-cross(w,Hvector);
    
    
    l_gain=1;
    mi=0.010;
    e0=0.01;freq=100;phi1=0;phi2=pi/2;phi3=pi;
    lamda=l_gain*exp(mi*det(A*A'));
    e1=e0*sin(freq*dur+phi1);
    e2=e0*sin(freq*dur+phi2);
    e3=e0*sin(freq*dur+phi3);
    E=abs(skew([e1 e2 e3]));
    A_hash=A'*inv(A*A'+lamda*E);
    d_dot=A_hash*Hdot;

    
    %%%%%%%%%%%% SATURATION %%%%%%%%%%%%%%
    ddmax=max(abs(d_dot));
    ddlim=deg2rad(50);
    ddgamma=ddmax/ddlim;
    if ddgamma>1  d_dot=d_dot/ddgamma; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % !!! BE CAREFUL FOR -torque !!!
    d=dt*d_dot+dprev;
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

    q_prev=q';
    w_prev=w;
    dd=[dd d_dot];
    TT=[TT Tc];
    
    if plot_motion==1 && mod(dur,35)==0
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
        object.vertices = [vnew(:,1) vnew(:,2) vnew(:,3)];
        object.faces = faces;
        patch(object,'FaceColor',[0.8 0.8 1.0],'EdgeColor','none','FaceLighting','gouraud','AmbientStrength', 0.15);
        camlight('headlight');
        material('dull');
        set(gca,'YTick',[]);set(gca,'XTick',[]);set(gca,'ZTick',[]);box on;grid on;
        axis(1*[-1 1 -1 1 -1 1]), hold off, view(-120,10)

        xlabel('x');ylabel('y');zlabel('z'); 
        pause(0.001); cla
    end

    
%     axes(SB2)
%     hold on
%     plot(dur,axiserror(3),'.red');
%     hold on
%     plot(dur,axiserror(2),'.green');
%     hold on
%     plot(dur,axiserror(1),'.blue');
%     if dur==endtime legend('error in axis z', 'error in axis y', 'error in axis x','Location','southeast') ,  end
%     axis([0 endtime -0.5 0.5]);
%     title('Error in x,y,z')
%     xlabel('time')
    
    
%     axes(SB3)
%     hold on
%     plot(dur,rad2deg(d1),'.red');
%     hold on
%     plot(dur,rad2deg(d2),'.green');
%     hold on
%     plot(dur,rad2deg(d3),'.blue');
%     hold on
%     plot(dur,rad2deg(d4),'.black');
%     if dur==endtime legend('d1', 'd2', 'd3','d4','Location','northeast') ,end
%     axis([0 endtime -270 270]);
%     title('Gimbal angles ')
%     xlabel('time')
%     ylabel('theta(degrees)')
%     yticks([-270:30:270]);
      d1all=[d1all d1];d2all=[d2all d2];d3all=[d3all d3];d4all=[d4all d4];
      axer=[axer ;axiserror(1) axiserror(2) axiserror(3)];
      manip(dur)=sqrt(det(A*A'));
      sigma=svd(A);
      lci_index(dur)=sigma(3)/sigma(1);
      hforplot(:,dur)=Hvector;
      wforplot(:,dur)=w;
end
%%
f1=figurem('Renderer','painters','Position',[0 0 1000 1600]);
xstyle='-black';ystyle='--black';zstyle=':black';

subplot(4,2,1)
hold on
plot(timeinsec,rad2deg(wforplot(1,:)),xstyle);
plot(timeinsec,rad2deg(wforplot(2,:)),ystyle);
plot(timeinsec,rad2deg(wforplot(3,:)),zstyle);
xticks([0:2:timeinsec(end)+1]);
yticks([-80:10:40])
yaxis(-80,40)
xlabel('Time, s')
ylabel('Angular Velocity, deg/s')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
t=title('a)','fontweight','normal')
set(t,'position',get(t,'position')-[11 150 0])
plot_dashed(singpoint,[-80:5:40],manip,dt)

% if approach==3
%     set(t,'position',get(t,'position')-[11 85 0])
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% elseif approach==4
%     set(t,'position',get(t,'position')-[11 85 0])
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% elseif approach==5 && approach2==2
%     set(t,'position',get(t,'position')-[11 85 0])
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% elseif approach==7 && approach2==2
%     set(t,'position',get(t,'position')-[11 125.5 0])
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% elseif approach==7 && approach2==3
%     set(t,'position',get(t,'position')-[11 118 0])
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% elseif approach==5 && approach2==3
%     set(t,'position',get(t,'position')-[11 125 0]) %+1/4
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% else
%     set(t,'position',get(t,'position')-[11 100 0]) %+1/4
%     plot_dashed(singpoint,[-80:5:40],manip,dt)
% end


subplot(4,2,2)
hold on
plot(timeinsec,hforplot(1,:),xstyle);
plot(timeinsec,hforplot(2,:),ystyle);
plot(timeinsec,hforplot(3,:),zstyle);
xticks([0:2:timeinsec(end)+1]);
xlabel('Time, s')
ylabel('CMG Momentum, Nms')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
plot_dashed(singpoint,[-1:0.1:3],manip,dt)
t=title('b)','fontweight','normal')
set(t,'position',get(t,'position')-[11 5 0])

subplot(4,2,3)
hold on
plot(timeinsec,rad2deg(axer(:,1)),xstyle);
plot(timeinsec,rad2deg(axer(:,2)),ystyle);
plot(timeinsec,rad2deg(axer(:,3)),zstyle);
ylabel('Quaternion Error, deg')
xlabel('Time, s')
legend('x', 'y','z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yticks([-10:10:50])
plot_dashed(singpoint,[-20:5:60],manip,dt)
t=title('c)','fontweight','normal')
set(t,'position',get(t,'position')-[11 97 0])

subplot(4,2,4)
hold on
plot(timeinsec,rad2deg(d1all),xstyle);
plot(timeinsec,rad2deg(d2all),ystyle);
plot(timeinsec,rad2deg(d3all),zstyle);
plot(timeinsec,rad2deg(d4all),'-.black');
ylabel('Gimbal Angles, deg')
xlabel('Time, s')
legend('d1', 'd2', 'd3','d4','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yticks([-180:45:180]);
plot_dashed(singpoint,[-190:190],manip,dt)
t=title('d)','fontweight','normal')
set(t,'position',get(t,'position')-[11 2*252 0])

% figure('Renderer','painters','Position',[20 50 1800 300])
subplot(4,2,5)
hold on
plot(timeinsec,rad2deg(dd(1,:)),xstyle)
plot(timeinsec,rad2deg(dd(2,:)),ystyle)
plot(timeinsec,rad2deg(dd(3,:)),zstyle)
plot(timeinsec,rad2deg(dd(4,:)),'-.black','linewidth',0.1)
ylabel('CMG Gimbal Rate, deg/s')
xlabel('Time, s')
legend('CMG_1', 'CMG_2', 'CMG_3','CMG_4','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
plot_dashed(singpoint,[-55:5:55],manip,dt)
yticks([-60:10:60]);
axis([0 20 -55 55])
t=title('e)','fontweight','normal')
set(t,'position',get(t,'position')-[11 140 0])


subplot(4,2,6)
hold on
plot(timeinsec,TT(1,:),xstyle)
plot(timeinsec,TT(2,:),ystyle)
plot(timeinsec,TT(3,:),zstyle)
ylabel('CMG Torque, Nm')
xlabel('Time, s')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
plot_dashed(singpoint,[-1.5:1:1.5],manip,dt)
axis([0 20 -1.5 1.5])
t=title('f)','fontweight','normal')
set(t,'position',get(t,'position')-[11 3.8 0])

subplot(4,2,7)
plot(timeinsec,manip,'-black')
ylabel('Manipulability Index')
xlabel('Time, s')
xticks([0:2:timeinsec(end)+1]);
plot_dashed(singpoint,[0:0.1:1.5],manip,dt)
t=title('g)','fontweight','normal')
set(t,'position',get(t,'position')-[11 1.8 0])

subplot(4,2,8)
plot(timeinsec,lci_index,'-black')
ylabel('Local Condition Index')
xlabel('Time, s')
xticks([0:2:timeinsec(end)+1]);
plot_dashed(singpoint,[0:0.1:1],manip,dt)
t=title('h)','fontweight','normal')
set(t,'position',get(t,'position')-[11 1.2 0])

[maxtorquex,ix]=max(abs(TT(1,:)));
[maxtorquey,iy]=max(abs(TT(2,:)));
[maxtorquez,iz]=max(abs(TT(3,:)));
maxtorques=[TT(1,ix) TT(2,iy) TT(3,iz)];
[mt,mti]=max(abs(maxtorques));
maxtorque=maxtorques(mti)
mti
disp('==========')
findsettling=0;
i=1;
while i<length(wforplot(:,stop+550:end)) && findsettling==0
    if abs(wforplot(1,stop+550+i))<deg2rad(0.001) && abs(wforplot(1,stop+550+i))<deg2rad(0.001) && abs(wforplot(1,stop+550+i))<deg2rad(0.001)
        findsettling=1;
        disp('settling time is:  ')
        stltime=(stop+550+i)/100 
        disp('sec')
    end
    i=i+1;
end
disp('=====min manip & singularity time===')
[min_manipulability,when]=min(manip)
disp('====min lci======')
disp('min lci');lci_min=min(lci_index)
disp('====manip settling======')
manip_settling=manip(end)
disp('=====LCI settling=====')
lci_settling=lci_index(end)
disp('======max wmega=========')
[maxwx,ix]=max(abs(wforplot(1,:)));
[maxwy,iy]=max(abs(wforplot(2,:)));
[maxwz,iz]=max(abs(wforplot(3,:)));
maxw=[wforplot(1,ix) wforplot(2,iy) wforplot(3,iz)];
[wmt,wmti]=max(abs(maxw));
maxw=maxw(wmti)
wmti
disp('======max momentum=========')
[maxhx,ix]=max(abs(hforplot(1,:)));
[maxhy,iy]=max(abs(hforplot(2,:)));
[maxhz,iz]=max(abs(hforplot(3,:)));
maxh=[hforplot(1,ix) hforplot(2,iy) hforplot(3,iz)];
[hmt,hmti]=max(abs(maxh));
maxh=maxh(hmti)
hmti


 MYTABLE=[MYTABLE;table(approach,approach2,approach3,maxtorque,mti,stltime,min_manipulability,when,lci_min,manip_settling,lci_settling,rad2deg(maxw),wmti,maxh,hmti)]
%   end
%%
fname = '/Users/charalamposp/Desktop/';
filename=strcat(num2str(approach),'_',num2str(approach2),'_',num2str(approach3))
saveas(gca, fullfile(fname, filename), 'eps');
approach
approach2
approach3
%   end
%%
