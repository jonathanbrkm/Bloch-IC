clear all; close all; clc;

om=1;
T=2*pi/om;

N0=1E3; 
N10=N0/2;
N20=N0/2;

beta10=pi; 
beta20=pi;

phi1=0*pi/180;
phi2=120*pi/180;

eta0_vec=linspace(0.9,1,200);
phi0_vec=linspace(0,2*pi,200);
[Eta0, Phi0]=meshgrid(eta0_vec,phi0_vec);

eta0_vec=Eta0(:);
phi0_vec=Phi0(:);

E00=1;
E0_vec=([0.8:0.01:0.9].*E00).';
delta=0;

eta_initial=[]; phi_initial=[];
eta_seg1=[]; phi_seg1=[];
eta_seg2=[]; phi_seg2=[];
eta_seg3=[]; phi_seg3=[];

for d=1:length(E0_vec)
d
    E0=E0_vec(d);
    
    theta=atan2(E0,delta); 

    N1=round(N10.*sqrt((E0.^2+delta.^2)./(E00.^2)));
    N2=round(N20.*sqrt((E0.^2+delta.^2)./(E00.^2)));

    beta1=beta10.*sqrt((E0.^2+delta.^2)./(E00.^2));
    beta2=beta20.*sqrt((E0.^2+delta.^2)./(E00.^2));

    N=N1+N2+N1;

    nx1=sin(theta)*cos(phi1);
    ny1=sin(theta)*sin(phi1);
    nz1=cos(theta);
    K1=[0,-nz1,ny1; nz1,0,-nx1; -ny1,nx1,0];
    Rx_180_60=expm(-(beta1/N1)*K1);

    nx2=sin(theta)*cos(phi2);
    ny2=sin(theta)*sin(phi2);
    nz2=cos(theta);
    K2=[0,-nz2,ny2;nz2,0,-nx2;-ny2,nx2,0];
    Ry_180_300=expm(-(beta2/N2)*K2);

    r(:,:,1)=[sqrt(1-eta0_vec.^2).*cos(phi0_vec),sqrt(1-eta0_vec.^2).*sin(phi0_vec),eta0_vec].';
    
    phi_initial=[phi_initial,atan2(r(2,:,1),r(1,:,1))];
    eta_initial=[eta_initial,r(3,:,1)];

    for f=1:N1
        r(:,:,f+1)=Rx_180_60*r(:,:,f);
    end
    eta_seg1=[eta_seg1,r(3,:,N1+1)];
    phi_seg1=[phi_seg1,atan2(r(2,:,N1+1),r(1,:,N1+1))];

    for f=(N1+1):(N1+N2)
        r(:,:,f+1)=Ry_180_300*r(:,:,f);
    end
    eta_seg2=[eta_seg2,r(3,:,N1+N2+1)];
    phi_seg2=[phi_seg2,atan2(r(2,:,N1+N2+1),r(1,:,N1+N2+1))];

    for f=(N1+N2+1):N
        r(:,:,f+1)=Rx_180_60*r(:,:,f);
    end
    eta_seg3=[eta_seg3,r(3,:,N+1)];
    phi_seg3=[phi_seg3,atan2(r(2,:,N+1),r(1,:,N+1))];

    clear r;
end

%%
phi_initial=mod(phi_initial,2*pi);
phi_seg1=mod(phi_seg1,2*pi);
phi_seg2=mod(phi_seg2,2*pi);
phi_seg3=mod(phi_seg3,2*pi);

Area_initial=(max(phi0_vec)-min(phi0_vec)).*(max(eta0_vec)-min(eta0_vec));
K_seg1=boundary(phi_seg1.',eta_seg1.');
Area_seg1=polyarea(phi_seg1(K_seg1),eta_seg1(K_seg1));
K_seg2=boundary(phi_seg2.',eta_seg2.');
Area_seg2=polyarea(phi_seg2(K_seg2),eta_seg2(K_seg2));
K_seg3=boundary(phi_seg3.',eta_seg3.');
Area_seg3=polyarea(phi_seg3(K_seg3),eta_seg3(K_seg3));

R32=Area_seg3./Area_seg2;
R21=Area_seg2./Area_seg1;
R10=Area_seg1./Area_initial;

R31=Area_seg3./Area_seg1;
R20=Area_seg2./Area_initial;

R30=Area_seg3./Area_initial;
%%
figure; set(gcf,'color','w'); 
plot(phi_initial,eta_initial,'.k',phi_seg1,eta_seg1,'.r',phi_seg2,eta_seg2,'.g',phi_seg3,eta_seg3,'.b');
grid on; axis([0 2*pi -1 1])
xlabel('$\phi$','interpreter','latex','fontsize',20);
ylabel('$\eta$','interpreter','latex','fontsize',20);
hold on;
plot(phi_initial,eta_initial,'.k',phi_seg1(K_seg1),eta_seg1(K_seg1),'--k',phi_seg2(K_seg2),eta_seg2(K_seg2),'--k',mod(phi_seg3(K_seg3),2*pi),eta_seg3(K_seg3),'--k');
legend('$0$','$1$','$2$','$3$','interpreter','latex','fontsize',14);
%%
figure;
sphere(25); view(150,20);
set(gcf,'color','w');
axis equal;
shading interp; colormap gray; alpha 0.25; 
line([-1 1],[0 0],[0 0],'LineWidth',1,'Color',[0 0 0]);
line([0 0],[-1 1],[0 0],'LineWidth',1,'Color',[0 0 0]);
line([0 0],[0 0],[-1 1],'LineWidth',1,'Color',[0 0 0]);
xlabel('$x$','interpreter','latex','fontsize',20,'position',[1.1 0 0]);
ylabel('$y$','interpreter','latex','fontsize',20,'position',[0 1.1 0]);
zlabel('$z$','interpreter','latex','fontsize',20,'position',[0 0 1.1]); h=gca; h.ZAxis.Label.Rotation=0;
hold on;

x_initial=sqrt(1-eta_initial.^2).*cos(phi_initial);
y_initial=sqrt(1-eta_initial.^2).*sin(phi_initial);
z_initial=eta_initial;

plot3(x_initial,y_initial,z_initial,'.k');
hold on;

x_seg1=sqrt(1-eta_seg1.^2).*cos(phi_seg1);
y_seg1=sqrt(1-eta_seg1.^2).*sin(phi_seg1);
z_seg1=eta_seg1;

plot3(x_seg1,y_seg1,z_seg1,'.r');
hold on;

x_seg2=sqrt(1-eta_seg2.^2).*cos(phi_seg2);
y_seg2=sqrt(1-eta_seg2.^2).*sin(phi_seg2);
z_seg2=eta_seg2;

plot3(x_seg2,y_seg2,z_seg2,'.g');
hold on;

x_seg3=sqrt(1-eta_seg3.^2).*cos(phi_seg3);
y_seg3=sqrt(1-eta_seg3.^2).*sin(phi_seg3);
z_seg3=eta_seg3;

plot3(x_seg3,y_seg3,z_seg3,'.b');
hold on;

%%
Etaf_tensor=reshape(eta_seg3,length(Phi0),length(Eta0),length(E0_vec));
Phif_tensor=reshape(phi_seg3,length(Phi0),length(Eta0),length(E0_vec));
Etai_tensor=reshape(eta_initial,length(Phi0),length(Eta0),length(E0_vec));
Phii_tensor=reshape(phi_initial,length(Phi0),length(Eta0),length(E0_vec));

G_tensor=NaN.*ones(length(Phi0),length(Eta0),length(E0_vec));

dphi0=2*pi./length(Phi0);
deta0=0.1./length(Phi0);
dE0=E0_vec(2)-E0_vec(1);

for i=2:length(E0_vec)-1
Etaf=Etaf_tensor(:,:,i); Etaf_prev=Etaf_tensor(:,:,i-1); Etaf_next=Etaf_tensor(:,:,i+1);
Phif=Phif_tensor(:,:,i); Phif_prev=Phif_tensor(:,:,i-1); Phif_next=Phif_tensor(:,:,i+1);
Etai=Etai_tensor(:,:,i);
Phii=Phii_tensor(:,:,i);

[dEtaf_dy,dEtaf_dx]=gradient(Etaf);
[dPhif_dy,dPhif_dx]=gradient(Phif);
detJf=(dPhif_dx.*dEtaf_dy)-(dPhif_dy.*dEtaf_dx);

[dEtai_dy,dEtai_dx]=gradient(Etai);
[dPhii_dy,dPhii_dx]=gradient(Phii);
detJi=(dPhii_dx.*dEtai_dy)-(dPhii_dy.*dEtai_dx);

dphif_dphii=((dPhif_dx.*dEtai_dy)-(dPhif_dy.*dEtai_dx))./detJi;
dphif_detai=((dPhif_dy.*dPhii_dx)-(dPhif_dx.*dPhii_dy))./detJi;
detaf_dphii=((dEtaf_dx.*dEtai_dy)-(dEtaf_dy.*dEtai_dx))./detJi;
detaf_detai=((dEtaf_dy.*dPhii_dx)-(dEtaf_dx.*dPhii_dy))./detJi;

detM=dphif_dphii.*detaf_detai-dphif_detai.*detaf_dphii;
dphif_dE0=(Phif_next-Phif_prev)./(2.*dE0);
detaf_dE0=(Etaf_next-Etaf_prev)./(2.*dE0);

detA=dphif_dphii.*detaf_dE0-dphif_dE0.*detaf_dphii;
detB=dphif_detai.*detaf_dE0-dphif_dE0.*detaf_detai;
G_mat=sqrt(detM.^2+detA.^2+detB.^2);

G_tensor(:,:,i)=G_mat;

% surf(G_mat);
% shading interp; colorbar; view(2);
% drawnow;
end

shpil=0;
G_tensor=G_tensor(1+shpil:end-shpil,1+shpil:end-shpil,:);
G_tensor=G_tensor(:);
G_avg=sum(G_tensor(~isnan(G_tensor))).*dphi0.*deta0.*dE0*(10/9);
G_norm=0.1*2*pi*0.1; 
G_avg_norm=G_avg./G_norm;

nexttile;
set(gcf,'color','w');
histogram(G_tensor(:),'FaceColor','g'); 
xlim([1 20]);
xlabel('$\mathcal{G}_{30}$','interpreter','latex','fontsize',20);
ylabel('$N(\mathcal{G}_{30})$','interpreter','latex','fontsize',20);
