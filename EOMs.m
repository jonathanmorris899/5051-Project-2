function [dX] = EOMs(t,X,vehicle_const,downburst_const)
    rm=downburst_const(1);
    zm=downburst_const(2);
    lambda=downburst_const(3);
    d=vehicle_const(1);
    mpayload=vehicle_const(2);
    rhogas=vehicle_const(3);
    Vol=4/3 * pi *(d/2)^3;
    Aref=pi*d.^2/4;
    g=9.81;
    x=X(1);
    y=X(2);
    z=X(3);
    vx=X(4);
    vy=X(5);
    vz=X(6);
    V=[vx,vy,vz];
    if z<0
        z=0;
        dX=[0,0,0,0,0,0]';
        return;
    end
    wind=calc_wind(x,y,z,rm,zm,lambda);
    [T, a, P, rho,nu,mu]=atmosisa(z);
    Fb=rho*g*Vol*[0,0,1];
    m=mpayload+rhogas*Vol;
    Fgrav=-m*g*[0,0,1];
    Air_rel_vel=V-wind;
    Re=rho*d*norm(Air_rel_vel)/(mu);
    Cd=calc_Cd(Re);
    drag_dir=Air_rel_vel/norm(Air_rel_vel);
    if norm(Air_rel_vel)==0
        Fdrag=0;
    else
        Fdrag=-Cd.*.5.*rho*(norm(Air_rel_vel)).^2.*Aref.*drag_dir;

    end
    a=(Fdrag+Fgrav+Fb)./m;
    % Fam=rho*g*Vol*a;
    % a=(Fdrag+Fgrav+Fb+Fam)./m;
    dX=[vx,vy,vz,a(1),a(2),a(3)]';

    % if t>6
    %     t
    %     Fdrag
    %     Fb
    %     Fgrav
    %     wind
    %     Air_rel_vel
    %     test=1;
    % end






end

function Cd=calc_Cd(Re)
    if Re==0
        Cd=0;
        return;
    end
    if Re<=1
        Cd=24/Re;
    elseif Re<=400 && Re>1
        Cd=24/(Re^(.646));

    elseif Re>400 && Re<(2.95*10^5)
        Cd=.5;
    elseif Re>2.95*10^5 && Re< 3.05*10^5
        Cd=.5+(.0803-.5)/(3.05*10^5-2.95*10^5)*(Re-2.95*10^5);
    elseif Re<=(2*10^6) && Re>(3.05*10^5)
        Cd=3.66*10^-4 *Re^(.4275);

    else
        Cd=.35;
        %error("invalid Re: " +num2str(Re))
    end

    


end

