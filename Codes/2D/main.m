clc;clear all,close all

format long

N=20;
[phi,del2_phi,del2_phi_c] = basis_functions(N);

A=[];
B=[];
C=[];
e=[];

parfor idx=1:184
    n = 2+ (idx-1)/10;
    A = A1_function (del2_phi,del2_phi_c, n, N);
    B = B1_function (phi,del2_phi,del2_phi_c, n, N);
    C = C1_function (phi, n, N);
    e = polyeig(A,B,C);
    e=real(e(imag(e) == 0));
    kappa=sort(sqrt(e));
    k0_N(idx) = kappa(1);

    A=[];
    B=[];
    C=[];
    e=[];
    kappa=[];
end
% % 
% k0 = 4.099276075446465;  % for n=1/10;
% k0 = 7.204007760688768;  % for n=1/3;
k0 = 4.159235819877813;  % for n=3;
% k0 = 1.849719369141224;  % for n=6;
% k0 = 1.296304915624706;  % for n=10;
% k0 = 1.166118732908618;    % for n=12;
% k0 = 0.881543660339067;    % for n=20;
% n=0.01:0.01:0.48;  
n=2:0.1:20.3;
plot(n,abs(k0*ones(size(n))-k0_N))
grid on
xlabel('refractive index n')
ylabel('error |k_0-k_0^N|')
title('Orginal n=3')
