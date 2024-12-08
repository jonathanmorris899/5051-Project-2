clear
clc


sz=1000;
Re=linspace(2.99*10^5,3.01*10^5,sz);

for i=1:length(Re)
    Cd(i)=calc_Cd(Re(i));
end

scatter(Re,Cd)
function Cd=calc_Cd(Re)
    if Re==0
        Cd=0;
        return;
    end
    if Re<=1
        Cd=24/Re;
    elseif Re<=400 && Re>1
        Cd=24/(Re^(.646));

    elseif Re>400 && Re<=(3*10^5)
        Cd=.5;
    elseif Re<=(2*10^6) && Re>(3*10^5)
        Cd=3.66*10^-4 *Re^(.4275);

    else
        Cd=.35;
    end




end