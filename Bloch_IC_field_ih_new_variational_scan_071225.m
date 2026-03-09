clear all; close all; clc;

om=1;
T=2*pi/om;

N0=1E2;

time_shift=linspace(-N0/10,N0/10,10); time_shift=0;
axis_shift=linspace(-pi/10,pi/10,10); axis_shift=0;

[T_mat,A_mat]=meshgrid(time_shift,axis_shift);

E00=1;
E0_vec=([0.8:0.01:0.9].*E00).';
delta=0;

for i=1:length(time_shift)
for j=1:length(axis_shift)
[i,j]
DT=T_mat(i,j);
DA=A_mat(i,j);

N10=N0/4+DT; N20=N0/2-2*DT;
tau10=T/4+(DT./N0)*T; tau20=T/2-(2*DT./N0)*T;
beta10=tau10; beta20=tau20;

phi1=0+DA; phi2=pi/2+DA; %co-rotating
% phi1=0+DA; phi2=pi/2-DA; %counter-rotating

eta0_vec=linspace(0.9,1,100);
phi0_vec=linspace(0,2*pi,100);
[Eta0,Phi0]=meshgrid(eta0_vec,phi0_vec);
eta0_vec=Eta0(:);
phi0_vec=Phi0(:);

eta_initial=[]; phi_initial=[];
eta_seg1=[]; phi_seg1=[];
eta_seg2=[]; phi_seg2=[];
eta_seg3=[]; phi_seg3=[];

for d=1:length(E0_vec)
E0=E0_vec(d);
theta=atan2(E0,delta);

N1=round(N10.*sqrt((E0.^2+delta.^2)./(E00.^2)));
N2=round(N20.*sqrt((E0.^2+delta.^2)./(E00.^2)));

beta1=beta10.*sqrt((E0.^2+delta.^2)./(E00.^2));
beta2=beta20.*sqrt((E0.^2+delta.^2)./(E00.^2));

N=N1+N2+N1;

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

R32(i,j)=Area_seg3./Area_seg2;
% R21=Area_seg2./Area_seg1;
% R10=Area_seg1./Area_initial;
% 
% R31=Area_seg3./Area_seg1;
% R20=Area_seg2./Area_initial;
% 
% R30=Area_seg3./Area_initial;

% figure; set(gcf,'color','w'); 
% plot(phi_initial,eta_initial,'.k',phi_seg1,eta_seg1,'.r',phi_seg2,eta_seg2,'.g',phi_seg3,eta_seg3,'.b');
% grid on; axis([0 2*pi -1 1])
% xlabel('$\phi$','interpreter','latex','fontsize',20);
% ylabel('$\eta$','interpreter','latex','fontsize',20);
% hold on;
% plot(phi_initial,eta_initial,'.k',phi_seg1(K_seg1),eta_seg1(K_seg1),'--k',phi_seg2(K_seg2),eta_seg2(K_seg2),'--k',mod(phi_seg3(K_seg3),2*pi),eta_seg3(K_seg3),'--k');
% legend('$0$','$1$','$2$','$3$','interpreter','latex','fontsize',14);
end
end

%%
figure; set(gcf,'color','w');
surf(time_shift./N0,axis_shift./pi,R32);
shading interp; colormap jet; view(2); colorbar;
xlabel('$\Delta\tau/T$','interpreter','latex','fontsize',20);
ylabel('$\Delta\phi/\pi$','interpreter','latex','fontsize',20);

