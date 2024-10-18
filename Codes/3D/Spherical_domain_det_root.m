clc;clear; close all;
format long
%% ------------------------------

  %Let D be a spherical domain of radius r 
  % with constant refractive index n(x) = n.

  % Hassan Yazdanin, 02/14/2023
%%
myfun = @(k,r,n,m) real(f(k,r,n,m)); % parameterized function

% parameter
r = 1;  %radius of the circule
n = 40;  % refractive index
m = 0;      %only m=0
x0 = 1.5; % initial point

fun = @(k) myfun(k,r,n,m);  % function of k alone
z = fzero(fun,x0)
