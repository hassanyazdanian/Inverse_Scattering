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
myfun = @(k,r,m) f_basis_eig (k,r,m); % parameterized function

% parameter
r = 1;  %radius of the circule
m = 9;
xl = 15; % lower
xu = 5;  %upper

fun = @(k) myfun(k,r,m);  % function of k alone
z = fzero(fun, xl)


% fplot(fun,[0,20])
% grid on

% m= 0;  z0,0= 3.196220616582541 z0,1=6.306437047688423  z0,2 = 9.439499137876405;  z0,3 = 12.577130640430653
% m= 1;  z1,0= 4.610899879049056 z1,1=7.799273800811232; z1,2 = 10.958067191919497  z1,3 = 14.108627805410926
% m= 2;  z2,0= 5.905678235420522 z2,1=9.196882599635321; z2,2 = 12.402220966864386
% m= 3;  z3,0= 7.143531023504839 z3,1=10.536669866589634; z3,2 = 13.795063594348296
% m= 4;  z4,0= 8.346605938750738 z4,1=11.836718456852012; z4,2 = 15.149870095505207
% m= 5;  z5,0= 9.525701355671723 z5,1=13.107363714984833
% m= 6;  z6,0= 10.687025855471223 z6,1=14.355156336629930
% m= 7;  z7,0= 11.834530212420384 z7,1=15.584551879279772
% m= 8;  z8,0= 12.970908653356057
% m= 8;  z9,0= 14.098093539271570