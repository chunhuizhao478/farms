% Fourth Order Runge Kutta Method to solve System of ODEs
% y1' = f1(ty1,y2)
% y2' = f2(ty1,y2)
% Fouth Order RK method for solving ODE
% y'= f(t.y)
% y_n+1 =y_n + (1/6) * ( delta _y1 + 2* delta_y2 + 2* delta _y3 + delta_y4)
% delta _y1 =h* f(tnyn)
% delta y2 = h * f(tn + h/2,yn + delta y1/2)
% delta_y3 =h * f(tn + h/2,yn + delta y2/2)
% delta_y4 = h * f(tn + hyn + delta_y3)
% dy1/dx = - 0.5 * y1
% dy2/dx = 4-0.3 * y2 - 0.1 * y1
% Numerically integrate above from x = 0 to x = 4 using step size of 0.5
% IC: y1(x=0) = 4; y2(×=0) = 6;
% viexact =
% vlexact =
clear
clc
clf
% Initial Values and time step
a = 0; % Lower limit of Integration
b = 2; % Upper limit of Integration
h= 0.05; % h = delta_x
y1(1) = 4;
y2(1) = 6;
n = (b-a)/h + 1; % no. of points
% Solution
x(1) = 0;
func1 = @ (y1) - 0.5 * y1; % function to integrate
func2 = @(v1, y2) 4 - 0.1 * y1 - 0.3 * y2; % function to integrate
% func_exact = @(x):
for i= 1:n-1
    x(i+1) = i*h;
    delta_y1_1 = h * feval(func1, y1(i));
    delta_y2_1 = h * feval(func1, y1(i) + (delta_y1_1)/2);
    delta_y3_1 = h * feval(func1, y1(i) + (delta_y2_1)/2);
    delta_y4_1 = h * feval(func1, y1(i) + delta_y3_1);
    y1(i+1) = y1(i) + (1/6) * ( delta_y1_1 + 2*delta_y2_1 + 2*delta_y3_1 + delta_y4_1);

    delta_y1_2 = h* feval(func2, y1(i), y2(1));
    delta_y2_2 = h* feval(func2, y1(i) + (delta_y1_1)/2, y2(i) + (delta_y1_2)/2);
    delta_y3_2 = h* feval(func2, y1(i) + (delta_y2_1)/2, y2(i) + (delta_y2_2)/2);
    delta_y4_2 = h* feval(func2, y1(i) + delta_y3_1    , y2(i) + delta_y3_2);
    y2(i+1) = y2(i) + (1/6) * ( delta_y1_2 + 2* delta_y2_2 + 2* delta_y3_2 + delta_y4_2);
end 
x
y1
y2
%{
yexact(1) = feval(func_exactx(1));
for i= 1n-1
x(i+1) = i*h;
yexact(i+1) = feval(func_exact,x(i+1));
end
yexact
%}
plot(x, y1, x, y2)
title("Plot of x vs y1 & y2 - System of ODE - 4th order Runge Kutta Method");
grid on
xlabel ("x")
ylabel ("y1, y2")
legend('y2','y2')
% Illustration for function with 2 variables
% al = 3
% b1 = 2
% fun1 = @(al,b1) al +