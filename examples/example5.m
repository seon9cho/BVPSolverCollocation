function example5

clc; close all;


options.algorithm = optimset('Display','off', 'Jacobian','off', ...
    'Algorithm','Levenberg-Marquardt','TolFun',1e-10);
    
% left boundary
a = 0;

% right boundary
b = 1;

% ODE
ode = @(x,y)(2*y);

% bc
bc = @(ya,yb) ya-1;

deg = 40
u0 = rand(deg+1,1);
    
[cf,f,exitflag] = fsolve(@(u) fun(u,a,b,deg,ode,bc),u0,options.algorithm);

dom = linspace(a,b,1000);
hold on;
plot(dom,exp(2*dom),'-k','LineWidth',2);
plot(dom,polyval(cf,dom),'--r','LineWidth',2);
h = xlabel('x');


% fsolve function
function out = fun(cf,a,b,deg,ode,bc)

    % nodes
    x = linspace(a,b,deg+3);

    % polynomial evaluated at nodes
    p = polyval(cf,x);
    
    % polynomial at the end points
    ya = p(1);
    yb = p(end);
    
    % coefficients of the derivative of the polynomial
    cf_der = cf(1:end-1).*fliplr(1:1:deg).';
    
    % deriatve of polynomial evaluated at the nodes
    p_der = polyval(cf_der,x);
        
    % output for fsolve
    out = [ bc(ya,yb), ... % boundary conditions
        p_der(2:end-1) - ode(x(2:end-1),p(2:end-1)) ... % ODE conditions
        ];

    






















