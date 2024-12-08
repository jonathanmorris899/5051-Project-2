clear
clc
close all

%%% 
%get into jet while vertical
%follow streamlines
%


umax=15;
D=2*10^3;
rm=1.1*D;
zm=.03*D;
lambda=umax/(.913*rm);

downburst_const=[rm,zm,lambda];


%% calc wind field more resolution
sz=50;
w=zeros(sz,sz,sz);
u=zeros(sz,sz,sz);
v=zeros(sz,sz,sz);
vmag=zeros(sz,sz,sz);
x_mat=linspace(-D*2.5,D*2.5,sz);
y_mat=linspace(-D*2.5,D*2.5,sz);
z_mat=linspace(0,1250,sz);
for i=1:sz
    for j=1:sz
        for k=1:sz
            x=x_mat(i);
            y=y_mat(j);
            z=z_mat(k);
            r=sqrt(x.^2+y.^2);

            wind=calc_wind(x,y,z,rm,zm,lambda);
            

            u(j,i,k)=wind(1);
            v(j,i,k)=wind(2);
            w(j,i,k)=-wind(3);
            vmag(j,i,k)=norm(wind);


        end
    end
end

w=-w;
[x,y,z]=meshgrid(x_mat,y_mat,z_mat);



%% simulation
X0=[0,-4000,900,0,0,0];

d=.625;
V=4/3*pi*(d/2)^3;

%mballoon=.128;
rho_gas=.166;






options=odeset('RelTol',1e-3,'AbsTol',1e-6);
tspan=[0,20*60];
[~,~,~,rho_air,~,~] = atmosisa(900);
mballoon=V*(rho_air-rho_gas)
vehicle_const=[d,mballoon,rho_gas];
[~,X]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);

[~,~,~,rho_air,~,~] = atmosisa(500);
mballoon=V*(rho_air-rho_gas)
vehicle_const=[d,mballoon,rho_gas];
[t2,X2]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);

[~,~,~,rho_air,~,~] = atmosisa(0);
mballoon=V*(rho_air-rho_gas)
vehicle_const=[d,mballoon,rho_gas];
[t3,X3]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);
%% slices 
figure()
hold on
slice(x_mat,y_mat,z_mat,w,[],[],[50,225,400])
colormap jet
c=colorbar();
c.Label.String="Wind Velocity [m/s]";
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
title("Vertical Wind Speed")
view(75,25)
% [xcyl,ycyl,zcyl]=cylinder(rm);
% surf(xcyl,ycyl,zcyl*500,'FaceAlpha',.10,'FaceColor','black')


figure()
slice(x_mat,y_mat,z_mat,vmag,[],[],[50,225,400])
colormap jet
c=colorbar();
c.Label.String="Wind Velocity [m/s]";
xlabel("X[m]")
ylabel("Y[m]")
zlabel("Z[m]")
title("Total Wind Speed")
view(75,25)



%% mass variation flow streamlines
N=12;
figure()
hold on
plot3(X(:,1),X(:,2),X(:,3),'LineWidth',2,'Color','k')
plot3(X2(:,1),X2(:,2),X2(:,3),'LineWidth',2,'Color','k','LineStyle','--')
plot3(X3(:,1),X3(:,2),X3(:,3),'LineWidth',2,'Color','k','LineStyle',':')

%potting downburst streamlines
% startx=cos(linspace(0,2*pi,N))*1300;
% starty=sin(linspace(0,2*pi,N))*1300;
% startz=ones(1,N)*400;
% streamline(x,y,z,u,v,w,startx,starty,startz,'Color','blue','HandleVisibility','off')
startx=cos(linspace(0,2*pi,N))*1500;
starty=sin(linspace(0,2*pi,N))*1500;
startz=ones(1,N)*400;
%streamline(x,y,z,u,v,w,startx,starty,startz,[.1,250],'Color','blue','HandleVisibility','off')
streamline(x,y,z,u,v,w,X(1,1),X(1,2),X(1,3),[.1,500],'Color','blue','LineWidth',2,'HandleVisibility','on')
startx=X2(500,1);
starty=X2(500,2);
startz=X2(500,3);
streamline(x,y,z,u,v,w,startx,starty,startz+13,[.1,500],'Color','blue','LineWidth',2,'HandleVisibility','on')
startx=X3(500,1);
starty=X3(500,2);
startz=X3(500,3);
streamline(x,y,z,u,v,w,startx,starty,startz,[.1,250],'Color','blue','LineWidth',2,'HandleVisibility','on')



xlim([-4500,4500])
ylim([-4500,4500])
zlim([0,1000])
[xcyl,ycyl,zcyl]=cylinder(rm);

%surf(xcyl,ycyl,zcyl*1,'FaceAlpha',1,'FaceColor','black')
L=legend("H_n = 900m, m=122g","H_n = 500m, m = 128g","H_n = 0m, m =135g","Wind Field Streamlines");


L.AutoUpdate='off';
vmag=sqrt(u.^2+v.^2+w.^2);

% t=slice(x_mat,y_mat,z_mat,vmag,[],[],60);
% colormap jet
% c=colorbar();
% c.Label.String="Wind Speed [m/s]";
% t.HandleVisibility='off';
view(90,0)
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
title("The Impact of Varying of Vehicle's Mass on Trajectory")




%% mass variation flow streamlines 2
figure()
hold on
startx=0;
starty=-4000;
startz=900;
plot3(X2(:,1),X2(:,2),X2(:,3),'LineWidth',2,'Color','k')
streamline(x,y,z,u,v,w,startx,starty,startz,[.1,10000],'Color','blue','LineWidth',2,'HandleVisibility','on')
startx=X2(400,1);
starty=X2(400,2);
startz=X2(400,3);
streamline(x,y,z,u,v,w,startx,starty,startz,[.1,300],'Color','blue','LineWidth',2,'HandleVisibility','off')
startx=X2(700,1);
starty=X2(700,2)+5;
startz=X2(700,3)+10;
streamline(x,y,z,u,v,w,startx,starty,startz,[.1,150],'Color','blue','LineWidth',2,'HandleVisibility','off')
xlim([-4500,4500])
ylim([-4500,4500])
zlim([0,1000])
view(90,0)
title("Trajectory of Balloon (H_n = 500 m) Compared to Flow Streamlines")
scatter3(X2(1,1),X2(1,2),X2(1,3),50,'red','diamond','filled')
plot3([0,0],[-5000,5000],[500,500],'--k')
legend("Vehicle Trajectory","Streamlines","Starting Location","H_n")


%% sensitivity to Location Y Hn=500m
%% sim
X0=[0,-9000,900,0,0,0];
[~,~,~,rho_air,~,~] = atmosisa(500);
mballoon=V*(rho_air-rho_gas)
vehicle_const=[d,mballoon,rho_gas];
[~,Xa]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines
X0=[0,-6000,900,0,0,0];
[~,Xb]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines
X0=[0,-3000,900,0,0,0];

[~,Xc]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines



%% plot
N=25;
startx=cos(linspace(0,2*pi,N))*1500;
starty=sin(linspace(0,2*pi,N))*1500;
startz=ones(1,N)*400;
figure()
hold on

plot3(Xa(:,1),Xa(:,2),Xa(:,3),'LineWidth',2,'Color','k')
plot3(Xb(:,1),Xb(:,2),Xb(:,3),'LineWidth',2,'Color','k','LineStyle','--')
plot3(Xc(:,1),Xc(:,2),Xc(:,3),'LineWidth',2,'Color','k','LineStyle',":")
scatter3(Xa(1,1),Xa(1,2),Xa(1,3),'red','diamond','filled')
scatter3(Xb(1,1),Xb(1,2),Xb(1,3),'red','diamond','filled','HandleVisibility','off')
scatter3(Xc(1,1),Xc(1,2),Xc(1,3),'red','diamond','filled','HandleVisibility','off')

streamline(x,y,z,u,v,w,Xa(697,1),Xa(697,2),Xa(697,3),[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','on')
streamline(x,y,z,u,v,w,Xb(555,1),Xb(555,2),Xb(555,3),[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','off')
streamline(x,y,z,u,v,w,Xc(415,1),Xc(415,2),Xc(415,3),[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','off')

[xcyl,ycyl,zcyl]=cylinder(rm);
surf(xcyl,ycyl,zcyl*500,'FaceAlpha',.10,'FaceColor','black')
legend("$Y_0$ = -9km","$Y_0$ = -6km","$Y_0$ = -3km","Starting Location","Wind Field Streamlines","Radius of W $\approx$ 0 m/s",'Interpreter','Latex')



view(45,45)
xlim([-10000,10000])
ylim([-10000,10000])
zlim([0,1000])
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
title("Varying $Y_0$ Trajectory Simulations",'Interpreter','Latex')


%% sensitivity to Location X Hn=500m
%% sim
X0=[0,-6000,900,0,0,0];
[~,~,~,rho_air,~,~] = atmosisa(500);
mballoon=V*(rho_air-rho_gas);
vehicle_const=[d,mballoon,rho_gas];
[~,Xa]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines
X0=[750,-6000,900,0,0,0];
[~,Xb]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines
X0=[1500,-6000,900,0,0,0];

[~,Xc]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);%potting downburst streamlines



%% plot
N=25;
startx=cos(linspace(0,2*pi,N))*1500;
starty=sin(linspace(0,2*pi,N))*1500;
startz=ones(1,N)*400;
figure()
hold on

plot3(Xa(:,1),Xa(:,2),Xa(:,3),'LineWidth',2,'Color','k')
plot3(Xb(:,1),Xb(:,2),Xb(:,3),'LineWidth',2,'Color','k','LineStyle','--')
plot3(Xc(:,1),Xc(:,2),Xc(:,3),'LineWidth',2,'Color','k','LineStyle',":")
scatter3(Xa(1,1),Xa(1,2),Xa(1,3),'red','diamond','filled')
scatter3(Xb(1,1),Xb(1,2),Xb(1,3),'red','diamond','filled','HandleVisibility','off')
scatter3(Xc(1,1),Xc(1,2),Xc(1,3),'red','diamond','filled','HandleVisibility','off')

streamline(x,y,z,u,v,w,Xa(555,1),Xa(555,2),Xa(555,3),[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','on')
streamline(x,y,z,u,v,w,Xb(587,1),Xb(587,2),Xb(587,3)+20,[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','off')
streamline(x,y,z,u,v,w,Xc(610,1),Xc(610,2),Xc(610,3),[.1,600],'Color','blue','LineWidth',2,'HandleVisibility','off')

[xcyl,ycyl,zcyl]=cylinder(rm);
surf(xcyl,ycyl,zcyl*500,'FaceAlpha',.10,'FaceColor','black')
legend("$X_0$ = 0km","$X_0$ = 0.75km","$X_0$ = 1.5km","Starting Location","Wind Field Streamlines","Radius of W $\approx$ 0 m/s",'Interpreter','Latex')



view(45,45)
xlim([-10000,10000])
ylim([-10000,10000])
zlim([0,1000])
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
title("Varying $X_0$ Trajectory Simulations",'Interpreter','Latex')







function [u,w]=velocity_field(r,z,rm,zm,lambda)
    c1=-.133;
    c2=1.1534;

    gamma=.85;
    delta=2;
    eps=2;
    kappa=0.6;
    chi=1.05;


    psi=(delta*r.^2/rm.^2).^gamma;



    u=lambda.*r./2 .*   (   exp(-(2.*gamma-psi).^2) + eps.*exp(-kappa.*(r.^2./rm.^2).^chi))  .* ( (z./zm) .^(c2-1) .* exp(c1*(z./zm).^c2));
    w=-lambda.* ( (1+2.*gamma.*psi.*(2.*gamma-psi)).*exp(-(2.*gamma-psi).^2) + eps.*exp(-kappa.*(r.^2./rm.^2).^chi) .* (1-kappa.*chi.*(r.^2/rm.^2).^chi))  .*...
        (zm./(c1.*c2).*(exp(c1.*(z./zm).^c2))-1);

end