clc;clear all;close all;
load('surface.mat');
MYTABLE=[];
%%%%%%%%%%%%%%%% SELECT APPROACH %%%%%%%%%%%%%%%%%%
approach=3
approach2=1
approach3=1

%  approach0=['1.0.0', '2.0.0', '3.0.0', '4.0.0', '5.1.0', '5.2.0', '5.3.0', '6.1.1', '6.2.1', '6.1.2', '6.2.2', '7.1.0', '7.2.0', '7.3.0', '7.4.0', '7.5.0'];
% for aaa=0:15
% approach=str2num(approach0(1+aaa*5))
% approach2=str2num(approach0(3+aaa*5))
% approach3=str2num(approach0(5+aaa*5))

% 1 - Generalized singularity robust inverse
% 2 - New method proposed by Bong Wie | 2005
% 3 - Standard pseudoinverse / Moore-Penrose
% 4 - Weighted pseudoinverse
% 5 - Preferred angles
%     5.1 with Generalized robust inverse
%     5.2 with Standard inverse
%     5.3 with Bong Wie|2005
% 6 - NULL space projection with the gradient of an index
%     6.1.1 Manipulability with GRI
%     6.2.1 LCI with GRI
%     6.1.2 Manipulability with BongWie
%     6.2.2 LCI with with BongWie
% 7 - Singularity robust inverse
%     7.1 Lamda=constant
%     7.2 Lamda with threshold | Manipulability
%     7.3 Lamda exponential | Manipulability
%     7.4 Lamda with threshold | LCI
%     7.5 Lamda exponential | LCI
% 8 - Conventional singularity with projection matrix
% 9 - Vadali/preferred angles (works only without Td)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=[0 0 0]';
J=1*[1 0 0; 0 1 0;0 0 1]; % Mass=24kg , Side of Cube=0.5m
Td=0*[0 10 0]';
Tc=[0 0 0]';
w_prev=w;
dt=0.01;
u=[0 0 0];
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
Kp=20/gfactor;
Ki=0.00001/gfactor;
Kwmega=100/gfactor;

gfactor=5;
Kp=2/gfactor;
Ki=0.00001/gfactor;
Kwmega=15/gfactor;
%%%%%%%%%%%%%%%%
endtime=30/dt;
timeinsec=[0:dt:endtime*dt-dt];
stop=3/dt;
plot_motion=0;

close all
%
SB1=figure('Renderer','painters','Position',[0 100 700 1600])
% SB2=subplot(1,3,2)
% SB3=subplot(1,3,3)
clc;
h0=1;
epsilon=deg2rad(1.5);
singpoint=0;
sum_delta_error=0;
Hdot=[0;0;0];
for dur=1:endtime
    %     if (d1<-pi/2+epsilon && d1>-pi/2-epsilon)&&(d3<pi/2+epsilon && d3>pi/2-epsilon) && (d2<0+epsilon && d2>0-epsilon) && (d4<0+epsilon && d4>0-epsilon)
    %         disp('singularity')
    %         singpoint=dur*dt;
    %     end
    
    if dur>stop
        Td=[0 0 0]';
        qref =[0.70711 -0.70711 0  0  ]; %
    end
    % expressed in the body frame
    wdot=inv(J)*(-skew(w)*J*w + Td +Tc); % the same as wdot=inv(J)*(-cross(w,(J*w)) + Td + Tc);
    w=wdot*dt+w_prev;
    
    qd=0.5*quatmultiply(q',[0;w]');
    if ( w(1)~=0 && w(2)~=0 && w(3)~=0 )
        q=quatmultiply(q_prev,[cos(norm(w)*dt/2);(w/norm(w))*sin(norm(w)*dt/2)]')';
    end
    att_angle(dur,:)=quat2eul(q');
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
    axer(dur,:)=u;
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
    
    if approach==1 % generalized singularity robust inverse
        l_gain=1;
        mi=0.10;
        e0=0.01;freq=100;phi1=0;phi2=pi/2;phi3=pi;
        lamda=l_gain*exp(-mi*det(A*A'));
        e1=e0*sin(freq*dur+phi1);
        e2=e0*sin(freq*dur+phi2);
        e3=e0*sin(freq*dur+phi3);
        E=abs(skew([e1 e2 e3])); % add +1*diag([1,1,1]) giati ta diagwnia stoixeia prepei na einai 1
        A_hash=A'*inv(A*A'+lamda*E);
        d_dot=A_hash*Hdot;
    elseif approach==2 % new method by bong wie | 2005
        l_gain=1e2;
        mi=1;
        e0=2e-1;
        freq=4e-4;
        phi1=pi/3;
        phi2=pi/2;
        phi3=pi;
        lamda=l_gain*exp(-mi*sqrt(det(A*A')));
        e1=e0*sin(freq*dur+phi1);
        e2=e0*sin(freq*dur+phi2);
        e3=e0*sin(freq*dur+phi3);
        V=lamda*abs(skew([e1 e2 e3])); % add +1*diag([1,1,1]) giati ta diagwnia stoixeia prepei na einai 1
        W=[1 lamda lamda lamda;
            lamda 2 lamda lamda;
            lamda lamda 3 lamda;
            lamda lamda lamda 4];
        A_hash=W*A'*inv(A*W*A'+V);
        d_dot=A_hash*Hdot;
    elseif approach==3 % Moore-Penrose
        A_hash=pinv(A);
        d_dot=A_hash*Hdot;
    elseif approach==4  % weighted pseudoinverse
        W=diag(0.01*[10 1 10 1]);
        A_hash=inv(W)*A'*inv(A*inv(W)*A');
        d_dot=A_hash*Hdot;
    elseif approach==5  % preferred angles
        if approach2==1  %preferred angles-generalized robust inverse
            gamma=0.1;
            l_gain=1;
            mi=0.010;
            e0=0.01;freq=100;phi1=0;phi2=pi/2;phi3=pi;
            lamda=l_gain*exp(mi*det(A*A'));
            e1=e0*sin(freq*dur+phi1);
            e2=e0*sin(freq*dur+phi2);
            e3=e0*sin(freq*dur+phi3);
            E=abs(skew([e1 e2 e3]));
            A_hash=A'*inv(A*A'+lamda*E);
        elseif approach2==2  %preferred angles-standard inverse
            gamma=2;
            A_hash=A'*inv(A*A');
        elseif approach2==3  %preferred angles-bong wie 2005
            gamma=0.1;
            l_gain=100;
            mi=0.010;
            e0=0.01;freq=100;phi1=0;phi2=pi/2;phi3=pi;
            lamda=l_gain*exp(mi*det(A*A'));
            e1=e0*sin(freq*dur+phi1);
            e2=e0*sin(freq*dur+phi2);
            e3=e0*sin(freq*dur+phi3);
            V=lamda*abs(skew([e1 e2 e3]));
            W=[1 lamda lamda lamda;
                lamda 2 lamda lamda;
                lamda lamda 3 lamda;
                lamda lamda lamda 4];
            A_hash=W*A'*inv(A*W*A'+V);
        end
        delta_desired=0*[pi/4 -pi/4 pi/4 -pi/4]';
        d_dot=A_hash*Hdot+gamma*(eye(4,4)-A_hash*A)*(delta_desired-[d1 d2 d3 d4]');
        %       gamma*(eye(4,4)-A_hash*A)*(delta_desired-[d1 d2 d3 d4]')
    elseif approach==6  % NULL space projection
        if approach2==1
            manip(dur)=sqrt(det(A*A'));
            NSgain=0.05;
            LCIFLAG=0;MANIPFLAG=1;
        elseif approach2==2
            sigma=svd(A);
            manip(dur)=sigma(3)/sigma(1);
            if approach3==1
                NSgain=0.8;
            elseif approach3==2
                NSgain=0.4;
            end
            LCIFLAG=1;MANIPFLAG=0;
        end
        deltastep=0.01;
        ALLCOMBIN=allcomb([d1 d1+deltastep d1-deltastep],[d2 d2+deltastep d2-deltastep],[d3 d3+deltastep d3-deltastep],[d4 d4+deltastep d4-deltastep]);
        for i=1:length(ALLCOMBIN)
            d1next=ALLCOMBIN(i,1);d2next=ALLCOMBIN(i,2);d3next=ALLCOMBIN(i,3);d4next=ALLCOMBIN(i,4);
            Anext=[-cos(b)*cos(d1next)   sin(d2next)             cos(b)*cos(d3next)   -sin(d4next);
                -sin(d1next)         -cos(b)*cos(d2next)     sin(d3next)          cos(b)*cos(d4next);
                sin(b)*cos(d1next)   sin(b)*cos(d2next)      sin(b)*cos(d3next)   sin(b)*cos(d4next)];
            sigma=svd(Anext);
            manext(i)=LCIFLAG*sigma(3)/sigma(1)+MANIPFLAG*sqrt(det(Anext*Anext'));
        end
        [maxw,maxw_index]=max(manext);
        delta_best=ALLCOMBIN(maxw_index,:);
        diafd=delta_best-[d1 d2 d3 d4];
        
        % %         exception handling
        % %         if manip(dur)>max(manext)
        % %             disp('wwwwwwpa');pause(2)
        % %         end
        
        for i=1:4 % NUMBER OF CMGS
            if diafd(i)~=0
                grad(i)=( maxw-manip(dur) )/diafd(i);
            else
                grad(i)=0;
            end
        end
        %%% the pinv gives bad/no results
        if approach3==1
            d_dot=GRI(A,dur)*Hdot+NSgain*(eye(4,4)-GRI(A,dur)*A)*grad';
        elseif approach3==2
            d_dot=bongwie(A,dur)*Hdot+NSgain*(eye(4,4)-bongwie(A,dur)*A)*grad';
        end
    elseif approach==7 % Singularity robust inverse
        if approach2==1
            lamda=0.1;
        elseif approach2==2
            l0=1;
            m=sqrt(det(A*A'));
            m0=0.1;
            makezero=1;
            if m>=m0 makezero=0; end
            lamda=makezero*l0*(1-m/m0)^2.0;
        elseif approach2==3
            l0=1;
            m=sqrt(det(A*A'));
            mi=0.2;
            lamda=l0*exp(-mi*m);
        elseif approach2==4
            l0=0.5;
            sigma=svd(A);
            m=sigma(3)/sigma(1);
            m0=0.1;
            makezero=1;
            if m>=m0 makezero=0; end
            lamda=makezero*l0*(1-m/m0)^2.0;
        elseif approach2==5
            l0=0.1;
            sigma=svd(A);
            m=sigma(3)/sigma(1);
            mi=0.2;
            lamda=l0*exp(-mi*m);
        end
        A_hash=A'*inv(A*A'+lamda*eye(3,3));
        d_dot=A_hash*Hdot;
    elseif approach==8 % Conventional singularity with projection matrix %%% EINAI VLAKEIA gt basika einai null space projection
        A_hash=A'*inv(A*A');
        d=[1 1 1 1]';
        p=0.5;
        n=(eye(4,4)-A_hash*A)*d;
        d_dot=A_hash*Hdot+p*n;
    elseif approach==9 % Vadali/Preferred angles
        delta_desired=[pi/4 -pi/4 pi/4 -pi/4]';
        a_gain=1;
        k_gain=1;
        d_dot=(eye(4,4)-A'*inv(A*A'+a_gain*eye(3,3))*A)*k_gain*delta_desired-[d1 d2 d3 d4]';
    end
    
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
    
    T1(:,dur)=cross(h(:,1),d1_vect);
    T2(:,dur)=cross(h(:,2),d2_vect);
    T3(:,dur)=cross(h(:,3),d3_vect);
    T4(:,dur)=cross(h(:,4),d4_vect);
    
    q_prev=q';
    w_prev=w;
    dd=[dd d_dot];
    TT=[TT Tc];
    
    if plot_motion==1 && mod(dur,15)==0
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
    eulerror(dur,:)=eulererror';
end
sum(abs(axer(:,2)))+sum(abs(axer(:,3)))

f1=figurem('Renderer','painters','Position',[0 0 1000 1600]);
xstyle='-black';ystyle='--black';zstyle=':black';text_X=-0.05;text_Y=-0.1;textsize=13; Nstep=300;

gca=subplot_tight(4,2,1)
hold on
plot(timeinsec,rad2deg(att_angle(:,3)),xstyle);
plot(timeinsec,rad2deg(att_angle(:,2)),ystyle);
plot(timeinsec,rad2deg(att_angle(:,1)),zstyle);
xlabel('Time, s')
ylabel('Attitude angle, deg')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=(yLimits(2)-yLimits(1))/Nstep;
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'a)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yticks([-100:10:-40, -30:20:50])

gca=subplot_tight(4,2,2)
hold on
plot(timeinsec,rad2deg(wforplot(1,:)),xstyle);
plot(timeinsec,rad2deg(wforplot(2,:)),ystyle);
plot(timeinsec,rad2deg(wforplot(3,:)),zstyle);
xlabel('Time, s')
ylabel('Angular Velocity, deg/s')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'b)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)


gca=subplot_tight(4,2,3)
hold on
plot(timeinsec,hforplot(1,:),xstyle);
plot(timeinsec,hforplot(2,:),ystyle);
plot(timeinsec,hforplot(3,:),zstyle);
xlabel('Time, s')
ylabel('CMG Momentum, Nms')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'c)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)



gca=subplot_tight(4,2,4)
hold on
% plot(timeinsec,rad2deg(axer(:,1)),xstyle);
% plot(timeinsec,rad2deg(axer(:,2)),ystyle);
% plot(timeinsec,rad2deg(axer(:,3)),zstyle);
plot(timeinsec,rad2deg(eulerror(:,3)),xstyle);
plot(timeinsec,rad2deg(eulerror(:,2)),ystyle);
plot(timeinsec,rad2deg(eulerror(:,1)),zstyle);
ylabel('Attitude Error, deg')
xlabel('Time, s')
legend('x', 'y','z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'d)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yaxis([-50 100])

gca=subplot_tight(4,2,5)
hold on
plot(timeinsec,rad2deg(d1all),xstyle);
plot(timeinsec,rad2deg(d2all),ystyle);
plot(timeinsec,rad2deg(d3all),zstyle);
plot(timeinsec,rad2deg(d4all),'-.black');
ylabel('Gimbal Angles, deg')
xlabel('Time, s')
legend('d1', 'd2', 'd3','d4','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'e)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)


% figure('Renderer','painters','Position',[20 50 1800 300])
gca=subplot_tight(4,2,6)
hold on
plot(timeinsec,rad2deg(dd(1,:)),xstyle)
plot(timeinsec,rad2deg(dd(2,:)),ystyle)
plot(timeinsec,rad2deg(dd(3,:)),zstyle)
plot(timeinsec,rad2deg(dd(4,:)),'-.black','linewidth',0.1)
ylabel('Gimbal Rates, deg/s')
xlabel('Time, s')
legend('CMG_1', 'CMG_2', 'CMG_3','CMG_4','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'f)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yaxis([-50 50])


gca=subplot_tight(4,2,7)
hold on
plot(timeinsec,TT(1,:),xstyle)
plot(timeinsec,TT(2,:),ystyle)
plot(timeinsec,TT(3,:),zstyle)
ylabel('CMG Torque, Nm')
xlabel('Time, s')
legend('x', 'y', 'z','Location','northoutside','orientation','horizontal')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'g)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
yaxis([-1 1])

gca=subplot_tight(4,2,8)
plot(timeinsec,manip,'-black')
ylabel('Manipulability Index')
xlabel('Time, s')
xticks([0:2:timeinsec(end)+1]);
yLimits = get(gca,'YLim');
limstep=((yLimits(2)-yLimits(1))/Nstep);
plot_dashed(singpoint,[yLimits(1):limstep:yLimits(2)],manip,dt);
text(text_X,text_Y,'h)','Units','Normalized','VerticalAlignment','Top','FontSize',textsize)
% yaxis([0 4])

% 
% %%
% subplot(4,2,8)
% plot(timeinsec,lci_index,'-black')
% ylabel('Local Condition Index')
% xlabel('Time, s')
% xticks([0:2:timeinsec(end)+1]);
% plot_dashed(singpoint,[0:0.1:1],manip,dt)
% t=title('h)','fontweight','normal')
% set(t,'position',get(t,'position')-[11 1.2 0])
% 
% [maxtorquex,ix]=max(abs(TT(1,:)));
% [maxtorquey,iy]=max(abs(TT(2,:)));
% [maxtorquez,iz]=max(abs(TT(3,:)));
% maxtorques=[TT(1,ix) TT(2,iy) TT(3,iz)];
% [mt,mti]=max(abs(maxtorques));
% maxtorque=maxtorques(mti)
% mti
% disp('==========')
% findsettling=0;
% i=1;
% while i<length(wforplot(:,stop+550:end)) && findsettling==0
%     if abs(wforplot(1,stop+550+i))<deg2rad(0.001) && abs(wforplot(1,stop+550+i))<deg2rad(0.001) && abs(wforplot(1,stop+550+i))<deg2rad(0.001)
%         findsettling=1;
%         disp('settling time is:  ')
%         stltime=(stop+550+i)/100
%         disp('sec')
%     end
%     i=i+1;
% end
% disp('=====min manip & singularity time===')
% [min_manipulability,when]=min(manip)
% disp('====min lci======')
% disp('min lci');lci_min=min(lci_index)
% disp('====manip settling======')
% manip_settling=manip(end)
% disp('=====LCI settling=====')
% lci_settling=lci_index(end)
% disp('======max wmega=========')
% [maxwx,ix]=max(abs(wforplot(1,:)));
% [maxwy,iy]=max(abs(wforplot(2,:)));
% [maxwz,iz]=max(abs(wforplot(3,:)));
% maxw=[wforplot(1,ix) wforplot(2,iy) wforplot(3,iz)];
% [wmt,wmti]=max(abs(maxw));
% maxw=maxw(wmti)
% wmti
% disp('======max momentum=========')
% [maxhx,ix]=max(abs(hforplot(1,:)));
% [maxhy,iy]=max(abs(hforplot(2,:)));
% [maxhz,iz]=max(abs(hforplot(3,:)));
% maxh=[hforplot(1,ix) hforplot(2,iy) hforplot(3,iz)];
% [hmt,hmti]=max(abs(maxh));
% maxh=maxh(hmti)
% hmti
% 
% 
% MYTABLE=[MYTABLE;table(approach,approach2,approach3,maxtorque,mti,stltime,min_manipulability,when,lci_min,manip_settling,lci_settling,rad2deg(maxw),wmti,maxh,hmti)]
% %   end

fname = '/Users/charalamposp/Desktop/';
filename=strcat(num2str(approach),'_',num2str(approach2),'_',num2str(approach3))
saveas(gca, fullfile(fname, filename), 'eps');
approach
approach2
approach3
%   end
%%
