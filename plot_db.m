clear
clc
close all

umax=20;
D=2*10^3;
rm=1.1*D;
zm=.03*D;
lambda=umax/(.913*rm);


sz=100;
w=zeros(sz,sz,sz);
u=zeros(sz,sz,sz);
v=zeros(sz,sz,sz);
%vh=zeros(sz,sz,sz);
x_mat=linspace(-D*2,D*2,sz);
y_mat=linspace(-D*2,D*2,sz);
z_mat=linspace(0,1000,sz);
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
            w(j,i,k)=wind(3);
            vmag(j,i,k)=norm(wind);


        end
    end
end
vrad=sqrt(u.^2+v.^2);
[x,y,z]=meshgrid(x_mat,y_mat,z_mat);
%% vert
figure()
slice(x_mat,y_mat,z_mat,w,[],[],[50,100,200,300,400])
colormap jet
c=colorbar();
c.Label.String="Wind Velocity";
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Vertical velocity slices downburst")


%% radial
figure()
slice(x_mat,y_mat,z_mat,vrad,[],[],[50,100,200,300,400])
colormap jet
c=colorbar();
c.Label.String="Wind Velocity";
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Radial Velocity slices downburst")

%% vector field
sz=25;
w=zeros(sz,sz,sz);
u=zeros(sz,sz,sz);
v=zeros(sz,sz,sz);
%vh=zeros(sz,sz,sz);
x_mat=linspace(-D*2,D*2,sz);
y_mat=linspace(-D*2,D*2,sz);
z_mat=linspace(0,1000,sz);
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
            w(j,i,k)=wind(3);
            vmag(j,i,k)=norm(wind);


        end
    end
end
vrad=sqrt(u.^2+v.^2);
[x,y,z]=meshgrid(x_mat,y_mat,z_mat);





figure ()
quiver3(x,y,z,u,v,w)
view(0,0)

%% 
figure()
