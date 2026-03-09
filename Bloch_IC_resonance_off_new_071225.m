clear all; close all; clc;

om=1;
T=2*pi/om;

N0=1E3;
N10=N0/4; N20=N0/2;
tau10=T/4; tau20=T/2;
beta10=tau10; beta20=tau20;
phi1=0; phi2=pi/2;

eta0_vec=linspace(0.9,1,400);
phi0_vec=linspace(0,2*pi,400);

[Eta0,Phi0]=meshgrid(eta0_vec,phi0_vec);
eta0_vec=Eta0(:);
phi0_vec=Phi0(:);

E0=1;
delta_vec=[0.4:0.01:0.6].';

eta_initial=[]; phi_initial=[];
eta_seg1=[]; phi_seg1=[];
eta_seg2=[]; phi_seg2=[];
eta_seg3=[]; phi_seg3=[];

for d=1:length(delta_vec)
    d
delta=delta_vec(d);
theta=atan2(E0,delta);

N1=round(N10.*sqrt((E0.^2+delta.^2)./(E0.^2)));
N2=round(N20.*sqrt((E0.^2+delta.^2)./(E0.^2)));

beta1=beta10.*sqrt((E0.^2+delta.^2)./(E0.^2));
beta2=beta20.*sqrt((E0.^2+delta.^2)./(E0.^2));
beta1_vec(d)=beta1;
beta2_vec(d)=beta2;

N=N1+N2+N1;
N_vec(d)=N;
N1_vec(d)=N1;
N2_vec(d)=N2;

nx1=sin(theta).*cos(phi1);
ny1=sin(theta).*sin(phi1);
nz1=cos(theta);
K1=[0,-nz1,ny1;nz1,0,-nx1;-ny1,nx1,0];

nx2=sin(theta).*cos(phi2);
ny2=sin(theta).*sin(phi2);
nz2=cos(theta);
K2=[0,-nz2,ny2;nz2,0,-nx2;-ny2,nx2,0];

Rx_90=expm(-(beta1./N1).*K1);
Ry_180=expm(-(beta2./N2).*K2);

r(:,:,1)=[sqrt(1-eta0_vec.^2).*cos(phi0_vec),sqrt(1-eta0_vec.^2).*sin(phi0_vec),eta0_vec].';
x_initial=r(1,:,1);
y_initial=r(2,:,1);
z_initial=r(3,:,1);

phi_initial=[phi_initial,atan2(y_initial,x_initial)];
eta_initial=[eta_initial,z_initial];

for f=1:N1-1
r(:,:,f+1)=Rx_90*r(:,:,f);
end

x_seg1=r(1,:,N1);
y_seg1=r(2,:,N1);
z_seg1=r(3,:,N1);

phi_seg1=[phi_seg1,atan2(y_seg1,x_seg1)];
eta_seg1=[eta_seg1,z_seg1];

for f=N1:(N1+N2)-1
r(:,:,f+1)=Ry_180*r(:,:,f);
end

x_seg2=r(1,:,N1+N2);
y_seg2=r(2,:,N1+N2);
z_seg2=r(3,:,N1+N2);

phi_seg2=[phi_seg2,atan2(y_seg2,x_seg2)];
eta_seg2=[eta_seg2,z_seg2];

for f=(N1+N2):(N1+N2+N1)-1
r(:,:,f+1)=Rx_90*r(:,:,f);
end

x_seg3=r(1,:,N);
y_seg3=r(2,:,N);
z_seg3=r(3,:,N);

phi_seg3=[phi_seg3,atan2(y_seg3,x_seg3)];
eta_seg3=[eta_seg3,z_seg3];


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
hold on;
P_3D=reshape(phi_seg3,sqrt(length(phi0_vec)),sqrt(length(phi0_vec)),length(delta_vec));
E_3D=reshape(eta_seg3,sqrt(length(phi0_vec)),sqrt(length(phi0_vec)),length(delta_vec));

set(gcf,'color','w'); box on; grid on; 
for k=1:length(delta_vec)
p_3d=P_3D(:,:,k); p_3d=p_3d(:);
e_3d=E_3D(:,:,k); e_3d=e_3d(:);
scatter3(p_3d,e_3d,delta_vec(k).*ones(length(p_3d),1),20,delta_vec(k).*ones(length(p_3d),1),'filled'); hold on;
end
alpha 0.5;
xlabel('$\phi$','interpreter','latex','fontsize',20);
ylabel('$\eta$','interpreter','latex','fontsize',20);
zlabel('$\Delta$','interpreter','latex','fontsize',20);
colorbar;

%%
Etaf_tensor=reshape(eta_seg1,length(Phi0),length(Eta0),length(delta_vec));
Phif_tensor=reshape(phi_seg1,length(Phi0),length(Eta0),length(delta_vec));
Etai_tensor=reshape(eta_initial,length(Phi0),length(Eta0),length(delta_vec));
Phii_tensor=reshape(phi_initial,length(Phi0),length(Eta0),length(delta_vec));

G_tensor=NaN.*ones(length(Phi0),length(Eta0),length(delta_vec));

dphi0=0.1./length(Phi0);
deta0=2*pi./length(Phi0);
ddelta=delta_vec(2)-delta_vec(1);

for i=2:length(delta_vec)-1
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
dphif_ddelta=(Phif_next-Phif_prev)./(2.*ddelta);
detaf_ddelta=(Etaf_next-Etaf_prev)./(2.*ddelta);

detA=dphif_dphii.*detaf_ddelta-dphif_ddelta.*detaf_dphii;
detB=dphif_detai.*detaf_ddelta-dphif_ddelta.*detaf_detai;
G_mat=sqrt(detM.^2+detA.^2+detB.^2);

G_tensor(:,:,i)=G_mat;

% surf(detM);
% shading interp; colorbar; view(2);
% drawnow;
end

shpil=3;
G_tensor=G_tensor(1+shpil:end-shpil,1+shpil:end-shpil,:);
G_tensor=G_tensor(:);
G_avg=sum(G_tensor(~isnan(G_tensor))).*dphi0.*deta0.*ddelta;
G_norm=0.1*2*pi*0.1;
G_avg_norm=G_avg./G_norm;

nexttile;
set(gcf,'color','w');
histogram(G_tensor(:),'FaceColor','g'); xlim([1 5]);
xlabel('$\mathcal{G}_{10}$','interpreter','latex','fontsize',20);
ylabel('$N(\mathcal{G}_{10})$','interpreter','latex','fontsize',20);

%%

figure;
set(gcf,'color','w');

yyaxis left;
h_left=plot(Rfi_vec,'o-');
h_left.Color='b';
h_left.MarkerFaceColor='b';
h_left.MarkerEdgeColor='k';

ax_left=gca;
ax_left.YColor='b'; 

ylabel('$R_{fi}$','interpreter','latex','fontsize',20,'color','b');

yyaxis right;
h_right=plot(Gfi_vec,'o-');
h_right.Color='r';
h_right.MarkerFaceColor='r';
h_right.MarkerEdgeColor='k'; 
ax_right=gca;
ax_right.YColor='r';
ylabel('$\bar{\mathcal G}_{fi}$','interpreter','latex','fontsize',20,'color','r');
xlabel('serial number','interpreter','latex','fontsize',20);


%%
phi=phi_initial;
eta=eta_initial;

Sigma_11=mean((phi-mean(phi)).^2);
Sigma_22=mean((eta-mean(eta)).^2);
Sigma_12=mean((phi-mean(phi)).*(eta-mean(eta)));
Sigma_21=Sigma_12;
Sigma_mat=[Sigma_11,Sigma_12;Sigma_21,Sigma_22];

lambda=eig(Sigma_mat);
det_cov=det(Sigma_mat);
trace_cov=Sigma_11+Sigma_22;
Area_cov=sqrt(lambda(1).*lambda(2));

%%
figure; set(gcf,'color','w'); box on;
hold on;

eta0_vec=linspace(0.9,1,400);
phi0_vec=linspace(0,2*pi,400);
[Eta0,Phi0]=meshgrid(eta0_vec,phi0_vec,delta_vec);

h=histogram2(Phi0,Eta0,'facecolor','k'); xlim([0 2*pi]); ylim([-1 1]);  shading interp;
set(h,'EdgeColor','none');

Q_phi=reshape(phi_seg1,length(Phi0),length(Eta0),length(delta_vec));
Q_eta=reshape(eta_seg1,length(Phi0),length(Eta0),length(delta_vec));
h=histogram2(Q_phi,Q_eta,'facecolor','r'); xlim([0 2*pi]); ylim([-1 1]); shading interp;
set(h,'EdgeColor','none');

Q_phi=reshape(phi_seg2,length(Phi0),length(Eta0),length(delta_vec));
Q_eta=reshape(eta_seg2,length(Phi0),length(Eta0),length(delta_vec));
h=histogram2(Q_phi,Q_eta,'facecolor','g'); xlim([0 2*pi]); ylim([-1 1]); shading interp;
set(h,'EdgeColor','none');

Q_phi=reshape(phi_seg3,length(Phi0),length(Eta0),length(delta_vec));
Q_eta=reshape(eta_seg3,length(Phi0),length(Eta0),length(delta_vec));
h=histogram2(Q_phi,Q_eta,'facecolor','b'); xlim([0 2*pi]); ylim([-1 1]); shading interp;
set(h,'EdgeColor','none');

xlabel('$\phi$','interpreter','latex','fontsize',20);
ylabel('$\eta$','interpreter','latex','fontsize',20);
view(-45,45);

%%
all_phi=[mod(phi_initial,2*pi); mod(phi_seg1,2*pi); mod(phi_seg2,2*pi); mod(phi_seg3,2*pi)];
all_eta=[eta_initial; eta_seg1; eta_seg2; eta_seg3];

[~,phi_grid,eta_grid]=histcounts2(all_phi,all_eta,'Normalization','pdf');

dphig=phi_grid(2)-phi_grid(1);
detag=eta_grid(2)-eta_grid(1);
dA=dphig.*detag;

rho_initial=histcounts2(mod(phi_initial,2*pi),eta_initial,phi_grid,eta_grid,'Normalization','pdf');
rho1=histcounts2(mod(phi_seg1,2*pi),eta_seg1,phi_grid,eta_grid,'Normalization','pdf');
rho2=histcounts2(mod(phi_seg2,2*pi),eta_seg2,phi_grid,eta_grid,'Normalization','pdf');
rho3=histcounts2(mod(phi_seg3,2*pi),eta_seg3,phi_grid,eta_grid,'Normalization','pdf');

p_init=sum(rho_initial.^2,'all')*dA;
p1=sum(rho1.^2,'all')*dA;
p2=sum(rho2.^2,'all')*dA;
p3=sum(rho3.^2,'all')*dA;

trace_rho_initial=p_init/p_init;
trace_rho1=p1/p_init;
trace_rho2=p2/p_init;
trace_rho3=p3/p_init;

%%
uncert_phi_initial=sqrt(mean(phi_initial.^2)-(mean(phi_initial)).^2);
uncert_eta_initial=sqrt(mean(eta_initial.^2)-(mean(eta_initial)).^2);
uncert_phi_seg1=sqrt(mean(phi_seg1.^2)-(mean(phi_seg1)).^2);
uncert_eta_seg1=sqrt(mean(eta_seg1.^2)-(mean(eta_seg1)).^2);
uncert_phi_seg2=sqrt(mean(phi_seg2.^2)-(mean(phi_seg2)).^2);
uncert_eta_seg2=sqrt(mean(eta_seg2.^2)-(mean(eta_seg2)).^2);
uncert_phi_seg3=sqrt(mean(phi_seg3.^2)-(mean(phi_seg3)).^2);
uncert_eta_seg3=sqrt(mean(eta_seg3.^2)-(mean(eta_seg3)).^2);

uncert_phi_vec=[uncert_phi_initial,uncert_phi_seg1,uncert_phi_seg2,uncert_phi_seg3];
uncert_eta_vec=[uncert_eta_initial,uncert_eta_seg1,uncert_eta_seg2,uncert_eta_seg3];

uncert_prod_initial=uncert_phi_initial.*uncert_eta_initial;
uncert_prod_seg1=uncert_phi_seg1.*uncert_eta_seg1;
uncert_prod_seg2=uncert_phi_seg2.*uncert_eta_seg2;
uncert_prod_seg3=uncert_phi_seg3.*uncert_eta_seg3;

uncert_prod_vec=[uncert_prod_initial,uncert_prod_seg1,uncert_prod_seg2,uncert_prod_seg3];

figure; set(gcf,'color','w'); box on;
hold on;
plot([0,1,2,3],uncert_phi_vec,'--ok','MarkerFaceColor',[0, 0.4470, 0.7410]);
plot([0,1,2,3],uncert_eta_vec,'--ok','MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot([0,1,2,3],uncert_prod_vec,'--ok','MarkerFaceColor',[0.4940, 0.1840, 0.5560]);
legend('$\sigma_{\phi}$','$\sigma_{\eta}$','$\sigma_{\phi}\sigma_{\eta}$','interpreter','latex','fontsize',20);
xlabel('Segment','interpreter','latex','fontsize',20);
xticks([0:1:3]);