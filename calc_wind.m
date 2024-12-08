function wind=calc_wind(x,y,z,rm,zm,lambda)
    if ~isreal(x) || ~isreal(y)
    disp("error")
    wind=[0,0,0]
    return;
    end
    theta=atan2(y,x);
    r=sqrt(x^2+y^2);


    c1=-.133;
    c2=1.1534;

    gamma=.85;
    delta=2;
    eps=2;
    kappa=0.6;
    chi=1.05;


    psi=(delta*r.^2/rm.^2).^gamma;



    vh=lambda.*r./2 .*   (   exp(-(2.*gamma-psi).^2) + eps.*exp(-kappa.*(r.^2./rm.^2).^chi))  .* ( (z./zm) .^(c2-1) .* exp(c1*(z./zm).^c2));
    w=lambda.* ( (1+2.*gamma.*psi.*(2.*gamma-psi)).*exp(-(2.*gamma-psi).^2) + eps.*exp(-kappa.*(r.^2./rm.^2).^chi) .* (1-kappa.*chi.*(r.^2/rm.^2).^chi))  .*...
        (zm./(c1.*c2).*(exp(c1.*(z./zm).^c2))-1);
    
    u=vh*cos(theta);
    v=vh*sin(theta);
    %u=u+25;
    
  
    wind=[u,v,w];
    if z>500
        wind=wind+[0,12.5,0];

    else
        wind=wind+[0,2.5,0];
    end
    
    %wind=[0,0,-10]
    
 


end