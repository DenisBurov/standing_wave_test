clear all
mu_0 = 4*pi*10^(-7);
a = 0.0525;
I=1;
x = linspace(-1, 1, 1000);
y = linspace(-1, 1, 1000);
z = linspace(-1, 1, 1000);
rho = sqrt(x.^2 + y.^2);
r = sqrt(x.^2 + y.^2 + z.^2);
alpha = sqrt(a.^2+r.^2-2.*a.*rho);
beta = sqrt(a.^2+r.^2+2.*a.*rho);
k = sqrt(1-alpha.^2./beta.^2);
gamma = sqrt(x.^2-y.^2);
C = mu_0.*I/pi;

% Legendre complete elliptic integrals of the first and second kind
m = 0.5; % parameter
d = a; %separation
[K,E] = ellipke(k.^2);
Bx = C.*x.*z/(2.*alpha.^2.*beta.*rho.^2).*((a.^2+r.^2).*E-alpha.^2.*K);
By = C.*y.*z/(2.*alpha.^2.*beta.*rho.^2).*((a.^2+r.^2).*E-alpha.^2.*K);
Bz = C./(2.*alpha.^2.*beta).*((a.^2-r.^2).*E+alpha.^2.*K);

% Plotting the magnetic field
figure;
plot(y, sqrt(Bx.^2+By.^2+Bz.^2)./max(sqrt(Bx.^2+By.^2+Bz.^2)), 'b', 'LineWidth', 2);
xlim([0 0.0525]);
xlabel('y, [m]');
ylabel('Normalized Magnetic Field Strength, [a.u.]');
title('Magnetic Field Strength vs. y');
grid on;


