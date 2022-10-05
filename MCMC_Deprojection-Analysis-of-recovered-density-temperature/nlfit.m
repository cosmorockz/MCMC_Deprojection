% this is a non-linear fit
% -------------------------------------------
clear all;
display('======================================');

kpc = 3.08e21;
keV = 1.16e7;
mu  = 0.598;
mue = 1.151;
kB  = 1.38e-16;

format shorte

model = @(b, x) b(1)./((x./b(2)).^b(3) + (x./b(2)).^b(4));

d = load('data.dat');
x = (d(:,1)+d(:,2))/(2.0);
ne = d(:,3);
T  = d(:,4)*keV;
prs = ne.*mue/mu*kB.*T;
beta0 = [1e-10 30.0 0.3 1.0];

[beta, R, J, Cov] = nlinfit(x, prs, model, beta0);
beta_err = sqrt(diag(Cov))';
r_prs = model(beta, x);

loglog(x, prs, 'o', x, r_prs)

% -------------------------------
% Saving the parameter data into text-file
% ------------------------------------
fileID = fopen('nlfit.param','w');
fprintf(fileID, '%4.4e ', beta);
fprintf(fileID, '\n');
fprintf(fileID, '%4.4e ', beta_err);
fclose(fileID);

A = [x r_prs];
fp = fopen('nlfit.out','w');
fprintf(fp, '%4.4e %4.4e\n', A');
fclose(fp);