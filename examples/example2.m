function exampl2


options.algorithm = optimset('Display','off', 'Jacobian','off', ...
    'Algorithm','Levenberg-Marquardt','TolFun',1e-10);
    

u0 = rand(2,1)
    
[u,f,exitflag] = fsolve(@(u) fun(u),u0,options.algorithm);

u

function out = fun(u)


out = [u(1)^2-4; 
       exp(u(2))-exp(-3);
       (u(1)-2)^2+(u(2)+3)^2+1e-8];

