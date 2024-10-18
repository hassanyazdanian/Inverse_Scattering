% function [phi,del2_phi,del2_phi_c] = basis_functions_3D_1(N)

clear all;close all;clc;
N=20;
syms r theta

z0_1 = 3.926602312047923;   %m = 0
z1_1 = 5.267657530336815;   %m = 1
z2_1 = 6.529929581924476;   %m = 2
z0_2 = 7.068582745628730;   %m = 0
z3_1 = 7.748589599213666;   %m = 3

z1_2 = 8.506949596505731;
z4_1 = 8.938692691525917;
z2_2 = 9.872617765818118; 
z5_1 = 10.108300420727264;
z0_3 = 10.210176122813033;

z3_2 = 11.190897323394417;
z6_1 = 11.262315022550924;
z1_3 = 11.687678377389791;
z7_1 = 12.403975429762559;
z4_2 = 12.475244230494308;

z2_3 = 13.104075687269328;
z0_4 = 13.351768777754097;
z8_1 = 13.535550975571837;
z5_2 = 13.733801310629357;
z3_3 = 14.476613461874880;

%Please see (4.1.9) and (4.1.10) in Nikolas Thesis
% \phi=0 and dphi/dr=0 on boundary

nc= -1.016170554350756e+04; %normalization_coeff

a1 = 35.874573920759836;
m=0;
phi(1) = (a1*spherical_jn(m, z0_1*r) + spherical_in(m, z0_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a1=-(-(2^(1/2)*pi^(1/2)*besseli(m + 1/2, z0_1*r)*(m + 1) - 2^(1/2)*z0_1*r*pi^(1/2)*besseli(m - 1/2, z0_1*r))/(2*z0_1*r^2*(1/(z0_1*r))^(1/2)))/(-(2^(1/2)*pi^(1/2)*besselj(m + 1/2, z0_1*r)*(m + 1) - 2^(1/2)*z0_1*r*pi^(1/2)*besselj(m - 1/2, z0_1*r))/(2*z0_1*r^2*(1/(z0_1*r))^(1/2)))
% a1=-(spherical_in(m-1, z0_1*r)-(m+1)/(z0_1*r)*spherical_in(m, z0_1*r))/(spherical_jn(m-1, z0_1*r)-(m+1)/(z0_1*r)*spherical_jn(m, z0_1*r))
%  a1=-d_spherical_in(m, z0_1*r)/d_spherical_jn(m, z0_1*r)

% del2_phi(1) = (-a1*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z0_1*r)*(2*m*sin(theta*(m - 1)) - z0_1^2*r^2*sin(theta*(m - 1)) + z0_1^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z0_1*r)*(2*m*sin(theta*(m - 1)) + z0_1^2*r^2*sin(theta*(m - 1)) - z0_1^2*r^2*sin(theta*(m + 1))))/(4*z0_1^2*r^4*sin(theta)*(1/(z0_1*r))^(3/2));
del2_phi(1) = del_2_phi_1(a1,m,z0_1);

% del2_phi_c(1)=conj(del2_phi(1));
del2_phi_c(1)=del2_phi(1);

a2 = 1.141240924219727e+02;
m=1;
phi(2) = (a2*spherical_jn(m, z1_1*r) + spherical_in(m, z1_1*r))*legendreP(m, cos(theta));
% a2=-d_spherical_in(m, z1_1*r)/d_spherical_jn(m, z1_1*r)

% del2_phi(2) = (-a2*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z1_1*r)*(2*m*sin(theta*(m - 1)) - z1_1^2*r^2*sin(theta*(m - 1)) + z1_1^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z1_1*r)*(2*m*sin(theta*(m - 1)) + z1_1^2*r^2*sin(theta*(m - 1)) - z1_1^2*r^2*sin(theta*(m + 1))))/(4*z1_1^2*r^4*sin(theta)*(1/(z1_1*r))^(3/2));
del2_phi(2) = del_2_phi_1(a2,m,z1_1);
del2_phi_c(2)=del2_phi(2);


a3=3.112685005085839e+02;
m=2; 
phi(3) = (a3*spherical_jn(m, z2_1*r) + spherical_in(m, z2_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a3=-d_spherical_in(m, z2_1*r)/d_spherical_jn(m, z2_1*r)

% del2_phi(3)=(-a3*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z2_1*r)*(2*m*sin(theta*(m - 1)) - z2_1^2*r^2*sin(theta*(m - 1)) + z2_1^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z2_1*r)*(2*m*sin(theta*(m - 1)) + z2_1^2*r^2*sin(theta*(m - 1)) - z2_1^2*r^2*sin(theta*(m + 1))))/(4*z2_1^2*r^4*sin(theta)*(1/(z2_1*r))^(3/2));
del2_phi(3) = del_2_phi_1(a3,m,z2_1);
del2_phi_c(3) = del2_phi(3);


a4=-8.304844085849090e+02;
m=0; 
phi(4) = (a4*spherical_jn(m, z0_2*r) + spherical_in(m, z0_2*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a4 = -d_spherical_in(m, z0_2*r)/d_spherical_jn(m, z0_2*r)

% del2_phi(4)=(-a4*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z0_2*r)*(2*m*sin(theta*(m - 1)) - z0_2^2*r^2*sin(theta*(m - 1)) + z0_2^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z0_2*r)*(2*m*sin(theta*(m - 1)) + z0_2^2*r^2*sin(theta*(m - 1)) - z0_2^2*r^2*sin(theta*(m + 1))))/(4*z0_2^2*r^4*sin(theta)*(1/(z0_2*r))^(3/2));
del2_phi(4) = del_2_phi_1(a4,m,z0_2);
del2_phi_c(4) = del2_phi(4);

a5=7.795940571831908e+02;
m=3;
phi(5) = (a5*spherical_jn(m, z3_1*r) + spherical_in(m, z3_1*r))*legendreP(m, cos(theta));

% dphi/dr=0 on boundary
% a5 = -d_spherical_in(m, z3_1*r)/d_spherical_jn(m, z3_1*r)

% del2_phi(5)=(-a5*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z3_1*r)*(2*m*sin(theta*(m - 1)) - z3_1^2*r^2*sin(theta*(m - 1)) + z3_1^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z3_1*r)*(2*m*sin(theta*(m - 1)) + z3_1^2*r^2*sin(theta*(m - 1)) - z3_1^2*r^2*sin(theta*(m + 1))))/(4*z3_1^2*r^4*sin(theta)*(1/(z3_1*r))^(3/2));
del2_phi(5) = del_2_phi_1(a5,m,z3_1);
del2_phi_c(5) = del2_phi(5);

a6 = -3.115411164836351e+03;
m = 1;
phi(6) = (a6*spherical_jn(m, z1_2*r) + spherical_in(m, z1_2*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a6 = -d_spherical_in(m, z1_2*r)/d_spherical_jn(m, z1_2*r)


% del2_phi(6)=(-a6*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z1_2*r)*(2*m*sin(theta*(m - 1)) - z1_2^2*r^2*sin(theta*(m - 1)) + z1_2^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z1_2*r)*(2*m*sin(theta*(m - 1)) + z1_2^2*r^2*sin(theta*(m - 1)) - z1_2^2*r^2*sin(theta*(m + 1))))/(4*z1_2^2*r^4*sin(theta)*(1/(z1_2*r))^(3/2));
del2_phi(6) = del_2_phi_1(a6,m,z1_2);
del2_phi_c(6) = del2_phi(6);

a7 = 1.847577188251642e+03;
m = 4;
phi(7) = (a7*spherical_jn(m, z4_1*r) + spherical_in(m, z4_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a7 = -d_spherical_in(m, z4_1*r)/d_spherical_jn(m, z4_1*r)

% del2_phi(7)=(-a7*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z4_1*r)*(2*m*sin(theta*(m - 1)) - z4_1^2*r^2*sin(theta*(m - 1)) + z4_1^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z4_1*r)*(2*m*sin(theta*(m - 1)) + z4_1^2*r^2*sin(theta*(m - 1)) - z4_1^2*r^2*sin(theta*(m + 1))))/(4*z4_1^2*r^4*sin(theta)*(1/(z4_1*r))^(3/2));
del2_phi(7) = del_2_phi_1(a7,m,z4_1);
del2_phi_c(7) = del2_phi(7);

a8 = -1.016170554350756e+04;
m = 2;
phi(8) = (a8*spherical_jn(m, z2_2*r) + spherical_in(m, z2_2*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a8 = -d_spherical_in(m, z2_2*r)/d_spherical_jn(m, z2_2*r)

% del2_phi(8)=(-a8*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z2_2*r)*(2*m*sin(theta*(m - 1)) - z2_2^2*r^2*sin(theta*(m - 1)) + z2_2^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z2_2*r)*(2*m*sin(theta*(m - 1)) + z2_2^2*r^2*sin(theta*(m - 1)) - z2_2^2*r^2*sin(theta*(m + 1))))/(4*z2_2^2*r^4*sin(theta)*(1/(z2_2*r))^(3/2));
del2_phi(8) = del_2_phi_1(a8,m,z2_2);
del2_phi_c(8) = del2_phi(8);
% 
a9 = 4.210383580919573e+03;
m = 5;
phi(9) = (a9*spherical_jn(m, z5_1*r) + spherical_in(m, z5_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a9 = -d_spherical_in(m, z5_1*r)/d_spherical_jn(m, z5_1*r)

% del2_phi(9)=(-a9*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z5_1*r)*(2*m*sin(theta*(m - 1)) - z5_1^2*r^2*sin(theta*(m - 1)) + z5_1^2*r^2*sin(theta*(m + 1))) ...
%                 -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z5_1*r)*(2*m*sin(theta*(m - 1)) + z5_1^2*r^2*sin(theta*(m - 1)) - z5_1^2*r^2*sin(theta*(m + 1))))/(4*z5_1^2*r^4*sin(theta)*(1/(z5_1*r))^(3/2));
del2_phi(9) = del_2_phi_1(a9,m,z5_1);
del2_phi_c(9) = del2_phi(9);

a10 =  1.921799834140692e+04;
m = 0;
phi(10) = (a10*spherical_jn(m, z0_3*r) + spherical_in(m, z0_3*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a10 = -d_spherical_in(m, z0_3*r)/d_spherical_jn(m, z0_3*r)

% del2_phi(10)=(-a10*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z0_3*r)*(2*m*sin(theta*(m - 1)) - z0_3^2*r^2*sin(theta*(m - 1)) + z0_3^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z0_3*r)*(2*m*sin(theta*(m - 1)) + z0_3^2*r^2*sin(theta*(m - 1)) - z0_3^2*r^2*sin(theta*(m + 1))))/(4*z0_3^2*r^4*sin(theta)*(1/(z0_3*r))^(3/2));
del2_phi(10) = del_2_phi_1(a10,m,z0_3);
del2_phi_c(10) = del2_phi(10);


a11 = -3.024822450897955e+04;
m = 3;
phi(11) = (a11*spherical_jn(m, z3_2*r) + spherical_in(m, z3_2*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a11 = -d_spherical_in(m, z3_2*r)/d_spherical_jn(m, z3_2*r)

% del2_phi(11)=(-a11*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z3_2*r)*(2*m*sin(theta*(m - 1)) - z3_2^2*r^2*sin(theta*(m - 1)) + z3_2^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z3_2*r)*(2*m*sin(theta*(m - 1)) + z3_2^2*r^2*sin(theta*(m - 1)) - z3_2^2*r^2*sin(theta*(m + 1))))/(4*z3_2^2*r^4*sin(theta)*(1/(z3_2*r))^(3/2));
del2_phi(11) = del_2_phi_1(a11,m,z3_2);
del2_phi_c(11) = del2_phi(11);

a12 = 9.316452337786426e+03;
m = 6;
phi(12) = (a12*spherical_jn(m, z6_1*r) + spherical_in(m, z6_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a12 = -d_spherical_in(m, z6_1*r)/d_spherical_jn(m, z6_1*r)

% del2_phi(12)=(-a12*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z6_1*r)*(2*m*sin(theta*(m - 1)) - z6_1^2*r^2*sin(theta*(m - 1)) + z6_1^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z6_1*r)*(2*m*sin(theta*(m - 1)) + z6_1^2*r^2*sin(theta*(m - 1)) - z6_1^2*r^2*sin(theta*(m + 1))))/(4*z6_1^2*r^4*sin(theta)*(1/(z6_1*r))^(3/2));
del2_phi(12) = del_2_phi_1(a12,m,z6_1);
del2_phi_c(12) = del2_phi(12);


a13 = 7.734408438586276e+04;
m = 1;
phi(13) = (a13*spherical_jn(m, z1_3*r) + spherical_in(m, z1_3*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a13 = -d_spherical_in(m, z1_3*r)/d_spherical_jn(m, z1_3*r)

% del2_phi(13)=(-a13*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z1_3*r)*(2*m*sin(theta*(m - 1)) - z1_3^2*r^2*sin(theta*(m - 1)) + z1_3^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z1_3*r)*(2*m*sin(theta*(m - 1)) + z1_3^2*r^2*sin(theta*(m - 1)) - z1_3^2*r^2*sin(theta*(m + 1))))/(4*z1_3^2*r^4*sin(theta)*(1/(z1_3*r))^(3/2));
del2_phi(13) = del_2_phi_1(a13,m,z1_3);
del2_phi_c(13) = del2_phi(13);

a14 = 2.014462473329677e+04;
m = 7;
phi(14) = (a14*spherical_jn(m, z7_1*r) + spherical_in(m, z7_1*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a14 = -d_spherical_in(m, z7_1*r)/d_spherical_jn(m, z7_1*r)

% del2_phi(14)=(-a14*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z7_1*r)*(2*m*sin(theta*(m - 1)) - z7_1^2*r^2*sin(theta*(m - 1)) + z7_1^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z7_1*r)*(2*m*sin(theta*(m - 1)) + z7_1^2*r^2*sin(theta*(m - 1)) - z7_1^2*r^2*sin(theta*(m + 1))))/(4*z7_1^2*r^4*sin(theta)*(1/(z7_1*r))^(3/2));
del2_phi(14) = del_2_phi_1(a14,m,z7_1);
del2_phi_c(14) = del2_phi(14);

a15 = -8.435488727102774e+04;
m = 4;
phi(15) = (a15*spherical_jn(m, z4_2*r) + spherical_in(m, z4_2*r))*legendreP(m, cos(theta));
% dphi/dr=0 on boundary
% a15 = -d_spherical_in(m, z4_2*r)/d_spherical_jn(m, z4_2*r)

% del2_phi(15)=(-a15*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z4_2*r)*(2*m*sin(theta*(m - 1)) - z4_2^2*r^2*sin(theta*(m - 1)) + z4_2^2*r^2*sin(theta*(m + 1))) ...
%                   -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z4_2*r)*(2*m*sin(theta*(m - 1)) + z4_2^2*r^2*sin(theta*(m - 1)) - z4_2^2*r^2*sin(theta*(m + 1))))/(4*z4_2^2*r^4*sin(theta)*(1/(z4_2*r))^(3/2));
del2_phi(15) = del_2_phi_1(a15,m,z4_2);
del2_phi_c(15) = del2_phi(15);

a16 = 2.765757652527587e+05;
m = 2;
phi(16) = (a16*spherical_jn(m, z2_3*r) + spherical_in(m, z2_3*r))*cos(m*theta);
% dphi/dr=0 on boundary
% a16 = -d_spherical_in(m, z2_3*r)/d_spherical_jn(m, z2_3*r)

del2_phi(16)=(-a16*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z2_3*r)*(2*m*sin(theta*(m - 1)) - z2_3^2*r^2*sin(theta*(m - 1)) + z2_3^2*r^2*sin(theta*(m + 1))) ...
                  -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z2_3*r)*(2*m*sin(theta*(m - 1)) + z2_3^2*r^2*sin(theta*(m - 1)) - z2_3^2*r^2*sin(theta*(m + 1))))/(4*z2_3^2*r^4*sin(theta)*(1/(z2_3*r))^(3/2));
del2_phi_c(16) = del2_phi(16);

a17 = -4.447177932366924e+05;
m = 0;
phi(17) = (a17*spherical_jn(m, z0_4*r) + spherical_in(m, z0_4*r))*cos(m*theta);
% dphi/dr=0 on boundary
% a17 = -d_spherical_in(m, z0_4*r)/d_spherical_jn(m, z0_4*r)

del2_phi(17)=(-a17*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z0_4*r)*(2*m*sin(theta*(m - 1)) - z0_4^2*r^2*sin(theta*(m - 1)) + z0_4^2*r^2*sin(theta*(m + 1))) ...
                  -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z0_4*r)*(2*m*sin(theta*(m - 1)) + z0_4^2*r^2*sin(theta*(m - 1)) - z0_4^2*r^2*sin(theta*(m + 1))))/(4*z0_4^2*r^4*sin(theta)*(1/(z0_4*r))^(3/2));
del2_phi_c(17) = del2_phi(17);

a18 = 4.275350605312533e+04;
m = 8;
phi(18) = (a18*spherical_jn(m, z8_1*r) + spherical_in(m, z8_1*r))*cos(m*theta);
% dphi/dr=0 on boundary
% a18 = -d_spherical_in(m, z8_1*r)/d_spherical_jn(m, z8_1*r)

del2_phi(18)=(-a18*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z8_1*r)*(2*m*sin(theta*(m - 1)) - z8_1^2*r^2*sin(theta*(m - 1)) + z8_1^2*r^2*sin(theta*(m + 1))) ...
                  -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z8_1*r)*(2*m*sin(theta*(m - 1)) + z8_1^2*r^2*sin(theta*(m - 1)) - z8_1^2*r^2*sin(theta*(m + 1))))/(4*z8_1^2*r^4*sin(theta)*(1/(z8_1*r))^(3/2));
del2_phi_c(18) = del2_phi(18);

a19 = -2.239463661041013e+05;
m = 5;
phi(19) = (a19*spherical_jn(m, z5_2*r) + spherical_in(m, z5_2*r))*cos(m*theta);
% dphi/dr=0 on boundary
% a19 = -d_spherical_in(m, z5_2*r)/d_spherical_jn(m, z5_2*r)

del2_phi(19)=(-a19*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z5_2*r)*(2*m*sin(theta*(m - 1)) - z5_2^2*r^2*sin(theta*(m - 1)) + z5_2^2*r^2*sin(theta*(m + 1))) ...
                  -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z5_2*r)*(2*m*sin(theta*(m - 1)) + z5_2^2*r^2*sin(theta*(m - 1)) - z5_2^2*r^2*sin(theta*(m + 1))))/(4*z5_2^2*r^4*sin(theta)*(1/(z5_2*r))^(3/2));
del2_phi_c(19) = del2_phi(19);

a20 = 9.084088013658640e+05;
m = 3;
phi(20) = (a20*spherical_jn(m, z3_3*r) + spherical_in(m, z3_3*r))*cos(m*theta);
% dphi/dr=0 on boundary
% a20 = -d_spherical_in(m, z3_3*r)/d_spherical_jn(m, z3_3*r)

del2_phi(20)=(-a20*2^(1/2)*pi^(1/2)*besselj(m + 1/2, z3_3*r)*(2*m*sin(theta*(m - 1)) - z3_3^2*r^2*sin(theta*(m - 1)) + z3_3^2*r^2*sin(theta*(m + 1))) ...
                  -2^(1/2)*pi^(1/2)*besseli(m + 1/2, z3_3*r)*(2*m*sin(theta*(m - 1)) + z3_3^2*r^2*sin(theta*(m - 1)) - z3_3^2*r^2*sin(theta*(m + 1))))/(4*z3_3^2*r^4*sin(theta)*(1/(z3_3*r))^(3/2));
del2_phi_c(20) = del2_phi(20);

%% compute Laplacian Numerically
% for i=1:N
%     del2_phi_1(i)=p_laplacian(phi(i));
%     del2_phi_c_1(i)=p_laplacian(conj(phi(i)));
% end

%% Plot basis functions
% FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\3D\Figures\Basis_Functions';   % Your destination folder
% 
% %  Compute derivative of basis functions
% for i=1:N
%     d_phi_r(i)=diff(phi(i),r);
% end
% 
% for i=1:N
%     figure(i);
%     r=0:0.01:1;
%     theta=241;
%     plot(r,eval(phi(i)),'LineWidth',1.5)
%     hold on
%     plot(r,eval(d_phi_r(i)),'r--','LineWidth',1.5)
%     grid on
%     leg1 = ['$\phi_{',num2str(i),'}$'];
%     leg2 = ['$d\phi_{',num2str(i),'}/dr$'];
%     legend(leg1,leg2,'Interpreter','Latex','FontSize', 15);
%     
%     ti=['$\theta = {',num2str(theta),'}$'];
%     title(ti,'Interpreter','Latex','FontSize', 15 );
%     xlabel('$r$','Interpreter','Latex','FontSize', 15 );
%     ylabel('$amplitude$','Interpreter','Latex','FontSize', 15 );
%     
%     FigName = ['phi_dpi_2d',num2str(i),'.fig'];
%     savefig(gcf,  fullfile(FolderName,FigName));
% 
%     FigName1 = ['phi_dpi_2d',num2str(i),'.emf'];
%     saveas(gcf,  fullfile(FolderName,FigName1))
% end
%% Comparison Laplacian Analytical vs Numerically
% for i=1:N
%     figure(i);
%     r=-1:0.01:1;
%     theta=241;
%     plot(eval(del2_phi(i)))
%     hold on
%     plot(eval(del2_phi_1(i)),'r--')
% end


% for i=1:N
%     figure(i);
%     r=-1:0.01:1;
%     theta=23;
%     plot(eval(del2_phi_c(i)))
%     hold on
%     plot(eval(del2_phi_c_1(i)),'r--')
% end
%% 3D plot of Basis functions
% FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\3D\Figures\Basis_Functions';   % Your destination folder
% % 
% for i=1:N
%     idx=0;
%     r=-1:0.01:1;
%     for theta=0:pi/100:2*pi
%         idx=idx+1;
%         phi_3d(:,idx)=eval(phi(i));
%     end
%     figure(i);
%     [Xi,Yi,Zi] = polarplot3d(phi_3d,'PlotType','surfc' ,'MeshScale' , [2 2]);
%     ti=['$\phi_{',num2str(i),'}$'];
%     title(ti,'Interpreter','Latex','FontSize', 15 );
%     xlabel('$x$','Interpreter','Latex','FontSize', 15 );
%     ylabel('$y$','Interpreter','Latex','FontSize', 15 );
%     zlabel('$z$','Interpreter','Latex','FontSize', 15 );
%     FigName = ['phi',num2str(i),'.fig'];
%     savefig(gcf,  fullfile(FolderName,FigName));
% 
%     FigName1 = ['phi',num2str(i),'.emf'];
%     saveas(gcf,  fullfile(FolderName,FigName1));
% end
%% 3D plot of deriviative of the Basis functions
FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\3D\Figures\Basis_Functions';   % Your destination folder
% 
%  Compute derivative of basis functions
for i=1:N
    d_phi_r(i)=diff(phi(i),r);
end

for i=1:N
    idx=0;
    r=-1:0.01:1;
    for theta=0:pi/100:2*pi
        idx=idx+1;
        d_phi_3d(:,idx)=eval(d_phi_r(i));
    end
    figure(i);
    [Xi,Yi,Zi] = polarplot3d(d_phi_3d,'PlotType','surfc' ,'MeshScale' , [2 2]);
    ti=['$d\phi_{',num2str(i),'}/dr$'];
    title(ti,'Interpreter','Latex','FontSize', 15 );
    xlabel('$x$','Interpreter','Latex','FontSize', 15 );
    ylabel('$y$','Interpreter','Latex','FontSize', 15 );
    zlabel('$z$','Interpreter','Latex','FontSize', 15 );
    FigName = ['dphi',num2str(i),'.fig'];
    savefig(gcf,  fullfile(FolderName,FigName));

    FigName1 = ['dphi',num2str(i),'.emf'];
    saveas(gcf,  fullfile(FolderName,FigName1));
end