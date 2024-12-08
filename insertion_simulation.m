clear
clc
close all
umax=15;
D=2*10^3;
rm=1.1*D;
zm=.03*D;
lambda=umax/(.913*rm);

downburst_const=[rm,zm,lambda];


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
szx=50;
szy=50;
w=-w;
[x,y,z]=meshgrid(x_mat,y_mat,z_mat);
xrect=linspace(-3000,-12000,szx);
yrect=linspace(-1250,1250,szy);
tspan=[0,20*60];
options=odeset('RelTol',1e-3,'AbsTol',1e-3  );
[~,~,~,rho_air,~,~] = atmosisa(500);
d=.625;
V=4/3*pi*(d/2)^3;
rho_gas=.166;
mballoon=V*(rho_air-rho_gas)
vehicle_const=[d,mballoon,rho_gas];
V=4/3*pi*(d/2)^3;

rho_gas=.166;
figure(1)
hold on
for i=1:szx
    i
    for j=1:szy
        
        
        X0=[yrect(j),xrect(i),900,0,0,0];
        X0_cell{i,j}=X0;
        [~,X_cell{i,j}]=ode45(@(t,X) EOMs(t,X,vehicle_const,downburst_const),tspan,X0,options);

    end
end
%%
close all
figure(1)
hold on
for i=1:szx
    i
    for j=1:szy
        
        
        X=X_cell{i,j};
        X0=X0_cell{i,j};
        %plot3(X(:,1),X(:,2),X(:,3))
        log(i,j)=check_inside(X,rm,250);
        if log(i,j) ==1
            
            scatter3(X0(1),X0(2),X0(3),'g','filled','HandleVisibility','off')
        else
            scatter3(X0(1),X0(2),X0(3),10,'r','filled','HandleVisibility','off')

        end
        

    end
end
%%
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")
[xcyl,ycyl,zcyl]=cylinder(rm);
zcyl=zcyl*250;
X0=X0_cell{25,25};
scatter3(X0(1),X0(2),X0(3),10,'g','filled')
X0=X0_cell{50,50};
scatter3(X0(1),X0(2),X0(3),10,'r','filled')
surf(xcyl,ycyl,zcyl,'FaceAlpha',.10,'FaceColor','black')
legend("Successful Insertion","Failed Insertion","Target Zone")
% plot3([-1250,-1250],[-3000,-12000],[900,900],'b')
% plot3([1250,1250],[-3000,-12000],[900,900],'b')
% 
% plot3([-1250,1250],[-3000,-3000],[900,900],'b')
% plot3([-1250,1250],[-12000,-12000],[900,900],'b')



xlabel("X[m]")
ylabel("Y[m]")
zlabel("Z[m]")

view(75,45)
ylim([-12500,5000])
xlim([-5000,5000])
zlim([0,1000])
title("Release Locations of Successful Downburst Insertions")





function log=check_inside(X,rcyl,hcyl)
    log=0;
    if min(X(:,3))> hcyl
        min(X(:,3));
        
        return;
    else
        for i=1:length(X)
            if X(i,1) <rcyl && X(i,1)>-rcyl
                if X(i,2) <rcyl && X(i,2)>-rcyl
                    log=1;
                    return;
                end

            end


        end
        


    end
    


end