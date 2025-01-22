clc; clear all; close all; beep off; 

ode = @(x,y) [y(2);-2*y(1).*(1-y(1).^2)];
jacobian = @(x,y)[0, 1; -2+6*y(1).^2, 0];

restpnt = [ 0,0; 
            1,0;
            -1,0];
      
time = 100;
scl = 3;
box = scl*[-1,1,-1,1];

plot_more = 'on';
num1 = 15;
num2 = 15;

phase_portrait(ode,jacobian,restpnt,box,time,plot_more,num1,num2)

x = linspace(-30,30,1000);
y = tanh(x);
z = 1-tanh(x).^2;

plot(y,z,'-g','LineWidth',2);

pt1 = [-1;0];
pt2 = [1;0];

dv1 = 0.5*[1;2];
dv2 = -0.5*[1;-2];

plot([pt1(1),pt1(1)+dv1(1)],[pt1(2),pt1(2)+dv1(2)],'--r','LineWidth',2);
plot([pt2(1),pt2(1)+dv2(1)],[pt2(2),pt2(2)+dv2(2)],'--r','LineWidth',2);



h = xlabel('u');
set(h,'FontSize',22);
h = ylabel('v');
set(h,'FontSize',22);








