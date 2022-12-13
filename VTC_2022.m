clc;clear all;close all;

img = imread('deck.png');     % Load a sample image
xImage = [-1.2 1.2; -1.2 1.2];  
yImage = [ 1.2 1.2; -1.2 -1.2];
zImage = [0 0; 0 0]; 


[vertices,faces,normals,name] = stlRead('Mini_Falcon_9_simple.stl');%satellite_keychain
cstep=1;
vx=vertices(1:cstep:end,1)-0;     %% prepei na paw to antikeimeno sto kentro
vy=vertices(1:cstep:end,2)-0;
vz=vertices(1:cstep:end,3)+abs(min(vertices(1:cstep:end,3)));
H=[vx;vy ;vz]; %Vertices of the cube
figure(1)
hold on
object.vertices = [vx vy vz]/141.732;
object.faces = [faces(1:cstep:end,1) faces(1:cstep:end,2) faces(1:cstep:end,3)];
patch(object,'FaceColor',       [0.8 0.8 1.0], ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
camlight('headlight');
material('dull');
view(20,10);axis equal;

filename='vtc2022.gif';
savetogif=0;
% vx=1*[-1 1 -1 1 -1 1 -1 1];
% vy=1*[-1 -1 1 1 -1 -1 1 1];
% vz=1*[-1 -1 -1 -1 1 1 1 1];
% H=[vx;vy ;vz]; %Vertices of the cube
% S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube
plot_motion=1;


q=[1 0 0 0];q_prev=q;
w=[0 0 0]';
J=diag(0.1*[1 1 1]);
M=diag([1 1 1]);
w_prev=[0;0;0];
v_prev=[0;0;0];
x_prev=[0;0;0];
T=[0;0;0];F=T;
dt=0.01;
qref=eul2quat([deg2rad(0) deg2rad(0) deg2rad(40)],'ZYX');
xref=[0;0;1];
k1=10;k2=10;
Integral=[0;0;0];
Kp=0.8;Ki=0.0001;Kwmega=0.8;
ur=[0;0;0];
stop=0;
endtime=30/dt;
timeinsec=[0:dt:endtime*dt-dt];
firsttime=1;
thx=0;thy=0;thz=0;
sumerr=[0;0;0];

ff=figure(1);
vv=0;
xy=[linspace(0,pi/2,endtime-800),linspace(pi/2,pi/2,800)];
thangle=[linspace(90,0,endtime/2), linspace(0,0,endtime/2)];


close all;ff=figure('Position',[0 100 800 800]);

for dur=1:endtime
    xref=[cos(xy(dur));cos(xy(dur));cos(xy(dur))];
    qref=eul2quat([deg2rad(0) deg2rad(0) deg2rad(thangle(dur))],'ZYX');
    W=[(rotz(thz)*roty(thy)*rotx(thx))'*[0;0;-0.1*9.81]]; % in the BF, to W dhmioyrgei roph alla egw pairnw to systhma to kentro mazas ara to cross(r,W)=0 opote den exw roph logw W
    acc=inv(M)*(F + W);
    wdot=inv(J)*(-skew(w)*J*w*0 + T);
    qd=0.5*quatmultiply(q,[0;w]');
    w=wdot*dt+w_prev;
    v=acc*dt+v_prev;
    if ( w(1)~=0 || w(2)~=0 || w(3)~=0 )
        q=quatmultiply(q_prev,[cos(norm(w)*dt/2);(w/norm(w))*sin(norm(w)*dt/2)]');
    end
    % Rotational error
    qerror=quatmultiply(quatconj(quatnormalize(qref)),quatnormalize((q)))';
    th=2*acos(qerror(1));
    axiserror=qerror(2:end);
    
    if qerror(1)<0
        axiserror=-1.0*axiserror;
    end
    
    Integral=Integral+axiserror;
    if dur==1
        axiserrorprev=axiserror;
    end
    axer(dur,:)=ur;
    ur=-Kp*axiserror-Ki*Integral-Kwmega*w;
    
    
    
    % Transl error
    x=v*dt+x_prev;
    
    xerr=(x-xref);
    sumerr=sumerr+xerr;
    if dur==1, xerrprev=xerr; end;
    ut=-8*xerr-8*(xerr-xerrprev)/dt-0.01*sumerr;
    xerrprev=xerr;
    
    
    FT=inv([1 0 0 0 1 ;
        0 1 0 1 0 ;
        0 0 1 0 0 ;
        0 k1 0 -k2 0 ;
        -k1 0 0 0 k2 ;
        ])*[ut;ur(1:2)];
    
    if dur<stop
        FT=[0;0;0;0;0];
    end
    
    
    if FT(3)<0, FT(3)=0; end
    
    F=[FT(1)+FT(5);FT(2)+FT(4);FT(3)];
    T=[k1*FT(2)-k2*FT(4);-k1*FT(1)+k2*FT(5);0];
    
    
    
    
    
    xmat(dur,:)=x;
    anglemat(dur,:)=quat2eul(q); thx=anglemat(dur,3);thy=anglemat(dur,2);thz=anglemat(dur,1);
    Tmat(dur,:)=T;
    Fmat(dur,:)=F;
    FTmat(dur,:)=FT;
    qerrmat(dur,:)=axiserror;
    
    
    x_prev=x;
    v_prev=v;
    q_prev=q;
    w_prev=w;
    
    
    VTCinert(:,dur)=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[FT(1) FT(2) FT(3)]';
    
    
    if plot_motion==1 && mod(dur,20)==0
        vnew=quatrotate(quatconj(q),[1 0 0]);
        quat=q;
        for i=1:length(vx)
            vnew(i,:)=quatrotate(quatconj(quat),[vx(i) vy(i) vz(i)]/141.732)+x'; % It uses the dcm matrix
        end
        
        
        
        subplot(2,3,[1 2 4 5])
        cla
        
        
        %         if firsttime==0
        %         cla(ff0);delete(ff0);
        %         delete(ft6); % arkei to delete na mpei sto telytaio gia na svhsei kai ta prohgoymena
        %         end
        %         firsttime=0;
        
        
        object.vertices = [vnew(:,1) vnew(:,2) vnew(:,3)];
        object.faces = faces;
        ff0=patch(object,'FaceColor',[0.6 0.6 1],'FaceLighting','gouraud','AmbientStrength', 0.15,'EdgeColor','none');
        camlight('headlight');
        material('dull');
        hold on
        
        f1=plot3(xmat(:,1),xmat(:,2),xmat(:,3),'--k')
        
        surf(xImage,yImage,zImage,'CData',img,'FaceColor','texturemap');

        axis([-1.2 1.2 -1.2 1.2 0 2]);
        %       axis([x(1)-0.2 x(1)+0.2 x(2)-0.2 x(2)+0.2 x(3)-0.2 x(3)+1.2]);
        %       axis equal;
        
        FT=-FT; % just to represent the direction of the expelled mass of each thruster
        
        lw=4;scalefactor=0.65;
        vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[FT(1) FT(2) FT(3)]';
        xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;0];
        ff1=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        ft1=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'  VTC');
        
        vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*W;
        xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;0.4];
        ff3=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        ft3=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'  Weight');
        
        scalefactor=2;
        vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0 FT(4) 0]';
        xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;1];
        ff2=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        ft2=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'  Thruster 1');
        
        vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[FT(5) 0 0]';
        xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;1];
        ff2=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        ft2=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'  Thruster 2');
        
        
        %
        %         vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[FT(4) 0 0]';
        %         xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;-1;0];
        %         ff4=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        %         ft4=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'Thruster 4');
        %
        %         vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0 FT(5) 0]';
        %         xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;1];
        %         ff5=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        %         ft5=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'Thruster 5');
        %
        %         vinert=scalefactor*(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0 FT(6) 0]';
        %         xn=(rotz(anglemat(dur,1))*roty(anglemat(dur,2))*rotx(anglemat(dur,3)))*[0;0;-1];
        %         ff6=quiver3(x(1)+xn(1),x(2)+xn(2),x(3)+xn(3),vinert(1),vinert(2),vinert(3),'linewidth',lw);
        %         ft6=text(x(1)+xn(1)+vinert(1),x(2)+xn(2)+vinert(2),x(3)+xn(3)+vinert(3),'Thruster 6');
        
        FT=-FT; % fix the sign
        
        pause(0.000000000000000001);
        xlabel('X');ylabel('Y');zlabel('Z');
        view(20+vv/5,10+vv/20);
        vv=vv+1;
        
        
        
        
        
        
        subplot(2,3,3);
        plot(timeinsec(1:dur),xmat);
        xlabel('Time, s');ylabel('Position, m'); legend('X','Y','Z','location','northeast');
        subplot(2,3,6);
        plot(timeinsec(1:dur),rad2deg(anglemat));
        xlabel('Time, s');ylabel('Attitude, deg'); legend('R','P','Y','location','northeast');
        %         subplot(2,3,5);
        %         plot(timeinsec(1:dur),Fmat);
        %         xlabel('Time, s');ylabel('Force, N'); legend('X','Y','Z');
        %         subplot(2,3,6);
        %         plot(timeinsec(1:dur),Tmat);ylabel('Torque, Nm'); legend('X','Y','Z');
        %         xlabel('Time, s');
        
        
%         if savetogif==1
%         frame = getframe(ff);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if firsttime==1
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%             firsttime=0;
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append');
%         end
%         end
        
        
    end
    
    
end

if savetogif==1
    for jj=0:250
        subplot(2,3,[1 2 4 5]);
        view(20+vv/5+jj/5,10+vv/20);
        frame = getframe(ff);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append');
        pause(1e-5);
    end
end



%%
plot3(xmat(:,1),xmat(:,2),xmat(:,3),'--k')

figurem,
subplot(2,3,1);
plot(timeinsec(1:dur),xmat);
xlabel('Time, s');ylabel('Position, m'); legend('X','Y','Z','location','southeast');
subplot(2,3,4);
plot(timeinsec(1:dur),rad2deg(anglemat));
xlabel('Time, s');ylabel('Attitude, deg'); legend('Y','P','R','location','southeast');
subplot(2,3,3);
plot(timeinsec(1:dur),Fmat);
xlabel('Time, s');ylabel('Force, N'); legend('X','Y','Z');
subplot(2,3,6);
plot(timeinsec(1:dur),Tmat);ylabel('Torque, Nm'); legend('X','Y','Z');


% theta=atan2(vecnorm(VTCinert(1:2,:)),VTCinert(3,:));
% phi=atan2(VTCinert(2,:),VTCinert(1,:));
%% gif player

gifplayer(filename,0.1)
[gifImage cmap] = imread(filename, 'Frames', 'all');
size(gifImage)
implay(gifImage);


