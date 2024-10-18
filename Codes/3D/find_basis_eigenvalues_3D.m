clc;clear; close all;
%% ------------------------------
 
format long

  %Let D be a spherical domain of radius r 
  % with constant refractive index n(x) = n.

  % Hassan Yazdanin, 02/14/2023
%%
myfun = @(k,r,m) real(f_basis_eig_3D (k,r,m)); % parameterized function
% myfun = @(k,r,m) f_basis_eig_3D (k,r,m); % parameterized function

% parameter
r = 1;  %radius of the circule
m = 0;
xl = 9; % lower
xu = 5;  %upper

fun = @(k) myfun(k,r,m);  % function of k alone
z = fzero(fun, xl)

% 
% fplot(fun,[0,9])
% grid on

% m= 0;  z0,0 = 3.926602312047923; z0,1 = 7.068582745628730; z0,2 = 10.210176122813033;  z0,3 = 13.351768777754097;
% m= 1;  z1,0 = 5.267657530336815; z1,1 = 8.506949596505731; z1,2 = 11.687678377389791;  z1,3 = 14.850490558473263;
% m= 2;  z2,0 = 6.529929581924476; z2,1 = 9.872617765818118; z2,2 = 13.104075687269328;
% m= 3;  z3,0 = 7.748589599213666; z3,1 = 11.190897323394417; z3,2 = 14.476613461874880;
% m= 4;  z4,0 = 8.938692691525917; z4,1 = 12.475244230494308; z4,2 = 15.815767135060669;
% m= 5;  z5,0 = 10.108300420727264; z5,1 = 13.733801310629357;
% m= 6;  z6,0 = 11.262315022550924; z6,1 = 14.971931193167560;
% m= 7;  z7,0 = 12.403975429762559; z7,1 = 16.193383255078270;
% m= 8;  z8,0 = 13.535550975571837;
% m= 9;  z9,0 = 14.658702887218029;