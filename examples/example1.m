function example1

clc; close all; clear all; beep off;

L = 20;

guess = @(x)[tanh(2*x);sech(x/2)];


options = bvpset('RelTol',1e-8,'AbsTol',1e-8);

solint = bvpinit(linspace(-L,0,30),guess);

sol = bvp5c(@ode,@bc,solint,options);

plot(sol.x,sol.y,'-k');

% -------------------------------------------------------------------------
% The ODE function
% -------------------------------------------------------------------------

function out = ode(x,y)

    u = y(1);
    v = y(2);
    
    out = [v;-2*u*(1-u^2)];

% -------------------------------------------------------------------------
% The boundary condition function
% -------------------------------------------------------------------------

function out = bc(ya,yb)

Proj = [-2,1];

out = [Proj*(ya-[-1;0]); yb(1)];











