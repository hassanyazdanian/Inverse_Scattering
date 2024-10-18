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
r = 1;  %radius of the circule
n = 6;  % refractive index

idx1 = 0;
for m=0:3             % Bessel function order
    idx1 = idx1+1;
    idx2 = 0;
    for k = 0:0.01:10  % k:eigenvalue
        idx2 = idx2+1;
        
        J1 = besselj(m,k*r); % Bessel function of first kind
        dJ1 = (m*besselj(m, k*r))/r - k*besselj(m + 1, k*r);
                
        J2  = besselj(m, k*n^(1/2)*r);
        dJ2 = (m*besselj(m, k*n^(1/2)*r))/r - k*n^(1/2)*besselj(m + 1, k*n^(1/2)*r);

        det(idx1,idx2) = J1* dJ2- J2*dJ1;
        
    end
end

k = 0:0.01:10;
plot(k, det(1,:),'b','LineWidth',2)

ti=['Distribution of eigenvalue for n = ',num2str(n),''];
title(ti,'Interpreter','Latex','FontSize', 15 );
xlabel(['wavenumber  ','$k$'],'Interpreter','Latex','FontSize', 15 );

grid on
% ylim([-1,1.1])
hold on
plot(k, det(2,:),'r-.','LineWidth',1.5)
plot(k, det(3,:),'k--','LineWidth',1.5)
plot(k, det(4,:),'g:','LineWidth',1.5)

legend('m=0','m=1','m=2','m=3','Interpreter','Latex','FontSize', 15,'Location','northwest')

FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\2D\Figures\Characterisitic_function';   % Your destination folder

FigName = ['Char_fun_n =',num2str(n),'.fig'];
savefig(gcf,  fullfile(FolderName,FigName));

FigName1 = ['Char_fun_n =',num2str(n),'.jpg'];
saveas(gcf,  fullfile(FolderName,FigName1));


% syms k m r n;
% J = besselj(m,k*r);
% g=diff(J, r);
% g

% syms k m r n;
% J = besselj(m,k*r);
% g=diff(J, r);
% g
% 
% j = sqrt(pi/(2*k*r))* besselj(m + 0.5, k*r);
% g1=diff(j, r)
% 
% y = sqrt(pi/(2* k*r))* besseli(m + 0.5, k*r);
% g2=diff(y, r)
