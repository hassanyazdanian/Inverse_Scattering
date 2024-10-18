clc;clear; close all;
%% ------------------------------
  % This codes is written based on part 5.1.1
  % of PhD thesis entitled "The Inverse Spectral 
  % Problem for the Reconstruction of the Refractive
  % Index from the Interior Transmission Problem" 
  % by Nikolaos Pallikarakis.


  %Let D be a circular domain of radius r 
  % with constant refractive index n(x) = n.

  % Hassan Yazdanin, 02/14/2023
%%
myfun = @(k,r,n,m) f(k,r,n,m); % parameterized function

% parameter
r = 1;  %radius of the circule
n = 10;  % refractive index
m = 0;
x0 = 0.5; % initial point

fun = @(k) myfun(k,r,n,m);  % function of k alone
z = fzero(fun,x0)
