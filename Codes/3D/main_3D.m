clc;clear all,close all

format long

N=8;


[phi,del2_phi,del2_phi_c] = basis_functions_3D_1(N);

% e = polyeig(A,-B,C);
A=[];
B=[];
C=[];
e=[];

% % % A = A_function (del2_phi,del2_phi_c, n, N);
% A = A1_function (del2_phi,del2_phi_c, n, N);
% % % B = B_function (phi,del2_phi,del2_phi_c, n, N);
% B = B1_function (phi,del2_phi,del2_phi_c, n, N);
% % % C = C_function (phi, n, N);
% C = C1_function (phi, n, N);
% 
% e = polyeig(A,B,C);
% e=real(e(imag(e) == 0));
% kappa=sort(sqrt(e))


% A1=A;A1(abs(A1)<1e-4)=0;
% B1=B;B1(abs(B1)<1e-4)=0;
% C1=C;C1(abs(C1)<1e-4)=0;
% e1 = polyeig(A1,B1,C1);
% e1=real(e1(imag(e1) == 0));
% kappa1=sort(sqrt(e1))
% 
% 
% A2=A;A2(abs(A2)<10)=0;
% B2=B;B2(abs(B2)<10)=0;
% C2=C;C2(abs(C2)<10)=0;
% e2 = polyeig(A2,B2,C2);
% e2=real(e2(imag(e2) == 0));
% kappa2=sort(sqrt(e2))
% n=30;
% idx=0;
% for n=0.1:0.01:0.5
% k0_N= zeros(14,1);
parfor idx=3:100
%     idx = idx+1;
% % n=50/100;
%     n = 1/idx;
%     n = 2+ (idx-1)/10;
%     n = 2*idx+1;
    n= idx;
    A = A1_function_3D (del2_phi,del2_phi_c, n, N);
    B = B1_function_3D (phi,del2_phi,del2_phi_c, n, N);
    C = C1_function_3D (phi, n, N);
    e = polyeig(A,B,C);
    e=real(e(imag(e) == 0));
    kappa=sort(sqrt(e));

%     A1=A;A1(abs(A1)<1e-4)=0;
%     B1=B;B1(abs(B1)<1e-4)=0;
%     C1=C;C1(abs(C1)<1e-4)=0;
%     e1 = polyeig(A1,B1,C1);
%     e1=real(e1(imag(e1) == 0));
%     kappa1=sort(sqrt(e1));

    k0_N(idx,1) = kappa(1);

    A=[];
    B=[];
    C=[];
    e=[];
    kappa=[];
end
%%
% k0 = 4.896642159292261;   % for n=1/6;
% n_true = 1/6;

% k0 = 5.070494255789265;  % for n=1/5;
% n_true = 1/5;

% k0 = 6.283180934909651;  % for n=1/4;
% n_true = 1/4;

% k0 = 7.696121917649135;    % for n=1/3;
% n_true = 1/3;

% k0 = 10.882999652679208; % for n=1/2;  !!!!
% n_true = 1/2;

% k0 = 7.695442854060307;  % for n=2;   !!!!
% n_true = 2;

% k0 = 4.443358060870907;  % for n=3;
% n_true = 3;

% k0 = 3.141596028192148;  % for n=4;
% n_true = 4;

% k0 = 2.267593967093401;  % for n=5;
% n_true = 5;

% k0 = 1.999045790544345;    % for n=6;
% n_true = 6;

% k0 = 1.816892308775226;    % for n=7;
% n_true = 7;

% k0 = 1.679747307928359;    % for n=8;
% n_true = 8;

% k0 = 1.570796326794897;    % for n=9;
% n_true = 9;

% k0 = 1.481207751924574;    % for n=10;
% n_true = 10;

% k0 = 1.405705333019371;    % for n=11;
% n_true = 11;

% k0 = 1.340874006175112;  % for n=12;
% n_true = 12;

% k0 = 1.190150476034732;  % for n=15;
% n_true = 15;

% k0 = 1.023434760320724;  % for n=20;
% n_true = 20;

% k0 = 0.830163979105108;  % for n=30;
% n_true = 30;

k0 = 0.716709033375646;  % for n=40;
n_true = 40;
 
% n=0.01:0.01:0.4;
n= 3:100;

% n=2*n+1;
plot(n,abs(k0*ones(size(n))-k0_N(3:100)')/k0,'LineWidth',1.5)
grid on
xlabel('refractive index $n$','Interpreter','Latex','FontSize', 15 );
ylabel('error $|k_0-k_0^N|$','Interpreter','Latex','FontSize', 15 );
ti=['Orginal $n=$ ', num2str(n_true)];
title(ti,'Interpreter','Latex','FontSize', 15 );
xline(n_true,':','LineWidth',2)
FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\3D\Figures\Estimation_n';   % Your destination folder
FigName = ['n=',num2str(n_true),'.fig'];
savefig(gcf,  fullfile(FolderName,FigName));
FigName1 = ['n=',num2str(n_true),'.emf'];
saveas(gcf,  fullfile(FolderName,FigName1))
