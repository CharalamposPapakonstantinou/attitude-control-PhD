
if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind);
end

%%
clear all ;clc;close all;
plotoption=0;

if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind);
end
n=[];c=[];qq=[];g=[];qu=[];er=[];gim=[];man=[];dd=[];g=[];
ini=[1 0 0];
vx=[-1 1 -1 1 -1 1 -1 1];
vy=[-1 -1 1 1 -1 -1 1 1];
vz=[-1 -1 -1 -1 1 1 1 1];
H=[vx;vy ;vz]; %Vertices of the cube
S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube

SB1=figurem;
k=1;
kq=1;
kg=1;
ke=1;
km=1;
kd=1;
kg=1;


s1 = serial('/dev/cu.usbmodem74508801','BaudRate',115200);
fopen(s1)
i=1;
plt=1;
while i<=10000
    r=fscanf(s1);
    pause(0.000001);
    
  

    
    if r(1)=='N'
        sn=strsplit(r(2:end),', ');
        n=[n ; str2num([sn{1},sn{2},sn{3}])];
        nt=deg2rad([n(end,3),n(end,2),n(end,1)]);
    elseif r(1)=='C'
        sc=strsplit(r(2:end),', ');
        c=[c ; str2num([sc{1},sc{2},sc{3}])];
        ct=deg2rad([c(end,3),c(end,2),c(end,1)]);
        
        if plotoption==0
%             subplot(4,2,6);
            plot(c(:,1),'-blue');
            plot(c(:,2),'-red');
            plot(c(:,3),'-green');
            hold on;
%             yaxis(-110,100);
            xlim([k-3,k+3]);
            k=k+1;
        else
            if mod(i,6)==0
                axes(SB1)
                q=eul2quat(ct,'ZYX')';
                quat=q';
                for i=1:length(vx)
                    vnew(i,:)=quatrotate(quatconj(quat),[vx(i) vy(i) vz(i)]); % It uses the dcm matrix
                    %               vnew(i,:)=(angle2dcm(ct(1),ct(2),ct(3))*[vx(i) vy(i) vz(i)]')';
                end
                H=vnew'; %Vertices of the cube
                S=[1 2 4 3; 1 2 6 5; 1 3 7 5; 3 4 8 7; 2 4 8 6; 5 6 8 7]; %Surfaces of the cube
                hold on
                
                set(gca,'visible','off')
                set(gca,'xtick',[],'ytick',[],'ztick',[]);
                for i=1:size(S,1)
                    Si=S(i,:);
                    ff=fill3(H(1,Si),H(2,Si),H(3,Si),[ i/size(S,1) 0 0],'facealpha',1);
                end
                
                axis([-3 3 -3 3 -3 3]), axis on, hold off;
                grid off;
                % xlabel('x');ylabel('y');zlabel('z');
                view(120,45);
                
            end
        end
        
 
    end
    
    
    
    
    i=i+1;
end

%%
m_min=min(g);
m_max=max(g);

m_min_x=m_min(1);
m_min_y=m_min(2);
m_min_z=m_min(3);
m_max_x=m_max(1);
m_max_y=m_max(2);
m_max_z=m_max(3);

biasx=(m_max_x+m_min_x)/2
biasy=(m_max_y+m_min_y)/2
biasz=(m_max_z+m_min_z)/2
%%
save('magnbiases.mat','biasx','biasy','biasz')

%%

close all;
g1=g(:,1);
g2=g(:,2);
g3=g(:,3);

hold on
scatter3(g1,g2,g3);
[c,r]=ellipsoid_fit([g1,g2,g3]);
% ellipsoid(c(1),c(2),c(3),r(2),r(1),r(3));
xlabel('x');ylabel('y');zlabel('z');
axis equal;
std([max(g1)-min(g1),max(g2)-min(g2),max(g3)-min(g3)])


[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( [g1,g2,g3] );
C = evecs*diag(1./radii)*evecs';
corp = (g - center')*C*2.4380e+03;  % to 2.64... prokyptei apo to 1/C(1,1)
corp1=corp(:,1);
corp2=corp(:,2);
corp3=corp(:,3);
scatter3(corp1,corp2,corp3);
std([max(corp1)-min(corp1),max(corp2)-min(corp2),max(corp3)-min(corp3)])

figurem,
subplot(1,3,1);hold on
scatter(g2,g3);
scatter(corp2,corp3);
axis equal
subplot(1,3,2);hold on
scatter(g1,g3);
scatter(corp1,corp3);
axis equal
subplot(1,3,3);hold on
scatter(g1,g2);
scatter(corp1,corp2);
axis equal
%%
% gc1=g1-(max(g1)+min(g1))/2;
% gc2=g2-(max(g2)+min(g2))/2;
% gc3=g3-(max(g3)+min(g3))/2;
%
% hold on
% scatter3(gc1,gc2,gc3);
% [cc,rr]=ellipsoid_fit([gc1,gc2,gc3]);
% % ellipsoid(cc(1),cc(2),cc(3),rr(2),rr(1),rr(3));
% axis equal
%%
% [A,b]=MinVolEllipse(g',0.01);
% gcal=(g-b')*A*8e6;
% factor=1;
% gcal1=gcal(:,1)*factor;
% gcal2=gcal(:,2)*factor;
% gcal3=gcal(:,3)*factor;
%
% hold on
% scatter3(gcal1,gcal2,gcal3);
% [ccal,rcal]=ellipsoid_fit([gcal1,gcal2,gcal3]);
% % ellipsoid(ccal(1),ccal(2),ccal(3),rcal(2),rcal(1),rcal(3));
% xlabel('x');ylabel('y');zlabel('z');
% axis equal;
% std([max(gcal1)-min(gcal1),max(gcal2)-min(gcal2),max(gcal3)-min(gcal3)])
%%

% [ofs,gain,rotM]=ellipsoid_fit_ST(g,0)
% [gain,rotM]=refine_3D_fit_ST(gain,rotM);
%
% gcal=(g-ofs')*rotM;
% factor=3000;
% gcal1=gcal(:,1)/(gain(1))*factor;
% gcal2=gcal(:,2)/(gain(2))*factor;
% gcal3=gcal(:,3)/(gain(3))*factor;
%
% hold on
% scatter3(gcal1,gcal2,gcal3);
% [ccal,rcal]=ellipsoid_fit([gcal1,gcal2,gcal3]);
% % ellipsoid(ccal(1),ccal(2),ccal(3),rcal(2),rcal(1),rcal(3));
% xlabel('x');ylabel('y');zlabel('z');
% axis equal;
%
% std([max(gcal1)-min(gcal1),max(gcal2)-min(gcal2),max(gcal3)-min(gcal3)])
%%
figurem,
subplot(1,3,1);hold on
scatter(g2,g3);
scatter(gc2,gc3);
scatter(gcal2,gcal3);
axis equal
subplot(1,3,2);hold on
scatter(g1,g3);
scatter(gc1,gc3);
scatter(gcal1,gcal3);
axis equal
subplot(1,3,3);hold on
scatter(g1,g2);
scatter(gc1,gc2);
scatter(gcal1,gcal2);
axis equal