clc;clear; close all;
%% ------------------------------
  


  %Let D be a spherical domain of radius r 
  % with constant refractive index n(x) = n.

  % Hassan Yazdanin, 09/17/2023
%%
r = 1;  %radius of the circule
n = 8;  % refractive index

idx1 = 0;
for m=0          %only m=0  % Bessel function order
    idx1 = idx1+1;
    idx2 = 0;
    for k = 0:0.01:10  % k:eigenvalue
        idx2 = idx2+1;
        
        j1 = spherical_jn(m,k*r); % Spherical Bessel function of first kind
        dj1 = pi^(1/2)*(k*besselj(m - 1/2, k*r) - (besselj(m + 1/2, k*r)*(m + 1/2))/r)*(1/(2*k*r))^(1/2) - (pi^(1/2)*besselj(m + 1/2, k*r))/(4*k*r^2*(1/(2*k*r))^(1/2));
%         dj1 =d_spherical_jn(m,k*r);
        
        j2  = spherical_jn(m, k*n^(1/2)*r);
        dj2 = - pi^(1/2)*((besselj(m + 1/2, k*n^(1/2)*r)*(m + 1/2))/r - k*n^(1/2)*besselj(m - 1/2, k*n^(1/2)*r))*(1/(2*k*n^(1/2)*r))^(1/2) - (pi^(1/2)*besselj(m + 1/2, k*n^(1/2)*r))/(4*k*n^(1/2)*r^2*(1/(2*k*n^(1/2)*r))^(1/2));
%         dj2 = d_spherical_jn(m, k*n^(1/2)*r);
        
        det(idx1,idx2) = j1* dj2- j2*dj1;
%       det(:,idx2) = (sqrt(n)*sin(k)*cos(sqrt(n)*k)-cos(k)*sin(sqrt(n)*k))/(sqrt(n)*k);
        
      %characteristic function for positive constant n not equal to 1 and
        %depend only on r (unit sphere located at center).
        % DOI: 10.1088/0266-5611/28/7/075005  eq. (17)
        % Theorem 3.1. Let η (sqrt(n)) be a constant not equal to either 0 or 1. The entire function F(k) (det) has all
        % real zeros when η is an integer or reciprocal of an integer. Otherwise, it has infinitely many
        % real and complex zeros.
    end
end
k = 0:0.01:10;
plot(k, det(1,:),'b','LineWidth',2)

ti=['Distribution of eigenvalue for n = ',num2str(n),''];
title(ti,'Interpreter','Latex','FontSize', 15 );
xlabel(['wavenumber  ','$k$'],'Interpreter','Latex','FontSize', 15 );

grid on
% ylim([-1,1.1])

legend('m=0','Interpreter','Latex','FontSize', 15,'Location','northeast')

FolderName = 'E:\University\PostDoc\Iran\Inverse_scattering\My_codes\3D\Figures\Characterisitic_function';   % Your destination folder

FigName = ['Char_fun_n =',num2str(n),'.fig'];
savefig(gcf,  fullfile(FolderName,FigName));

FigName1 = ['Char_fun_n =',num2str(n),'.emf'];
saveas(gcf,  fullfile(FolderName,FigName1));


% syms k m r n;
% J = besselj(m,k*r);
% g=diff(J, r);
% g

% syms k m r n;
% i = sp_besseli(m,k*r);
% g=diff(i, r);
% g