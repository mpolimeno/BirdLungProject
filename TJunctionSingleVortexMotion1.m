function TJunctionSingleVortexMotion1(cf)

k = sqrt(5);

mean_minus = 0;
amp_minus = 1;
mean_plus = 0.1;
amp_plus = 0.75;

rtol = 1e-6;
atol = 1e-9;
options = odeset('RelTol',rtol,'AbsTol',atol);

%%% HARD CODE initial vortex position in zeta plane %%%
% load('TJunctionPotential_3-2.mat')
% 
% zeta_vort = -1/2 + 1i + 0.05 + 0.02*1i;
% z_Kutta = -1; % corner at which to enforce Kutta condition (+1 or -1)
% tspan = [0 pi];
% 
% zeta_vort = 1/2 + 1i + 0.05 - 0.02*1i;
% z_Kutta = 1;
% tspan = [0 pi];
% 
% zeta_vort = 1/2 + 1i - 0.05 - 0.02*1i;
% z_Kutta = 1;
% tspan = [pi 2*pi];
% 
% zeta_vort = -1/2 + 1i - 0.05 - 0.02*1i;
% z_Kutta = -1;
% tspan = [pi 2*pi];
% 
% if (imag(zeta_vort) > 1) 
%     ReZeta_vort = ReZeta2;
%     ImZeta_vort = ImZeta2;
%     Z_vort = Z2;
% else
%     ReZeta_vort = ReZeta1;
%     ImZeta_vort = ImZeta1;
%     Z_vort = Z1;
% end;
% 
% x_vort = interp2(ReZeta_vort,ImZeta_vort,real(Z_vort),real(zeta_vort),imag(zeta_vort));
% y_vort = interp2(ReZeta_vort,ImZeta_vort,imag(Z_vort),real(zeta_vort),imag(zeta_vort));
% z_IC = x_vort + 1i*y_vort;

z_Kutta = 1;
dt = 10^(-2);
%tspan = [dt pi];
tspan = [pi+dt 2*pi];
cyc_dir = sign(sin(tspan(1))); % 1 in forward cycle; -1 in reverse; was cos previously, but should give same answer

if (mean_minus == 0 && mean_plus == 0)
    kappa_lead = pi/2*(1 - 1/k^2)*cyc_dir*(amp_minus*(z_Kutta - 1/k) - amp_plus*(z_Kutta + 1/k));
    z_IC = z_Kutta + dt*sign(kappa_lead)*abs(kappa_lead)^(1/2)/(4*sqrt(2))*10^(1/4)*(sqrt(5/3) + sign(kappa_lead)*1i); 
else
    kappa_lead = pi/2*(1 - 1/k^2)*(mean_minus*(z_Kutta - 1/k) - mean_plus*(z_Kutta + 1/k));
    z_IC = z_Kutta + sqrt(dt)*sign(kappa_lead)*abs(kappa_lead)^(1/2)/4*10^(1/4)*(sqrt(5/3) + sign(kappa_lead)*1i);
end;
zeta_IC = TJunctionPotential_Zeta_Func(z_IC,k);

disp(['vortex z-position is z = ' num2str(z_IC) '; zeta = ' num2str(zeta_IC)])
%pause % FOUND IT!

[T,Z] = ode45(@(T,Z) RHSFunc(T,Z,k,mean_minus,amp_minus,mean_plus,amp_plus,z_Kutta),tspan,z_IC,options);

zeta_vec = NaN*ones(size(Z));
for ind = 1:length(T)
    zeta_vec(ind) = TJunctionPotential_Zeta_Func(Z(ind),k);
end;

figure(cf)
subplot(4,2,[1 3])
hold on

% construct lines showing T-junction boundary
line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
line([.5 2],[1 1],'LineWidth',3,'Color','k')
line([-2 2],[0 0],'LineWidth',3,'Color','k')

plot(real(zeta_vec),imag(zeta_vec),'.-b','MarkerSize',16)

hold off

xlabel('Re($\zeta$)','interpreter','latex')
ylabel('Im($\zeta$)','interpreter','latex')

subplot(4,2,2)
plot(T,real(zeta_vec),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($\zeta$)','interpreter','latex')

subplot(4,2,4)
plot(T,imag(zeta_vec),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($\zeta$)','interpreter','latex')

subplot(4,2,[5 7]); 
plot(real(Z),imag(Z),'.-b','MarkerSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')

subplot(4,2,6)
plot(T,real(Z),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($z$)','interpreter','latex')

subplot(4,2,8)
plot(T,imag(Z),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($z$)','interpreter','latex')

function out = RHSFunc(t,z,k,mean_minus,amp_minus,mean_plus,amp_plus,z_Kutta)
% right-hand side of vortex motion equation: dz/dt = ...

u_minus = u_func(t,mean_minus,amp_minus);
u_plus = u_func(t,mean_plus,amp_plus);

gamma_vort = abs(z_Kutta - z)^2/imag(z)*(-u_minus/(z_Kutta + 1/k) + u_plus/(z_Kutta - 1/k));
%gamma_vort = 0;

%rhs = 1/abs(Fprime(z,k))^2*(1/pi*(u_minus/(z + 1/k) - u_plus/(z - 1/k)) + 1i*gamma_vort/(2*pi*(z - conj(z))) + 1i*gamma_vort/(4*pi)*F2prime(z,k)/Fprime(z,k));
rhs = 1/abs(Fprime(z,k))^2*(1/pi*(u_minus/(z + 1/k) - u_plus/(z - 1/k)) + 1i*gamma_vort/(2*pi*(z - conj(z))) + 1i*gamma_vort/(4*pi)*F2prime_over_Fprime(z,k));

out = conj(rhs);

function out = u_func(t,a,b)

out = a + b*sin(t);
% out = b; % CHANGE BACK LATER

function out = Fprime(z,k)
% derivative of conformal map from upper half-plane to T-junction

out = 5/pi*sqrt(1 - z.^2)./(1 - (k*z).^2);

function out = F2prime(z,k)
% second derivative of conformal map from upper half-plane to T-junction

out = -5/pi*z./sqrt(1 - z.^2).*(1 + k^2*(z.^2 - 2))./(1 - (k*z).^2).^2;

function out = F2prime_over_Fprime(z,k)

out = -z/(1 - z^2)*(1 + k^2*(z^2 - 2))/(1 - (k*z)^2);