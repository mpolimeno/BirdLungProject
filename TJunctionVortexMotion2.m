function [T,zeta_mat] = TJunctionVortexMotion2(tpe,cf)
% tpe = 1 for decoupled equations for each vortex
% tpe = 2 for coupled equations for vortex velocity & Kutta condition

k = sqrt(5);

mean_minus = 0;
amp_minus = 1*1;
mean_plus = 0.;
amp_plus = 1*0.75;

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

z_IC = NaN*ones(2,1);
zeta_IC = NaN*ones(2,1);
dt = 10^(-2);
%tspan = [dt pi];
tspan = [pi+dt 2*pi];
cyc_dir = sign(sin(tspan(1))); % 1 in forward cycle; -1 in reverse

for ind = 1:2 % loop over two corners
    z_Kutta = 2*ind - 3; % ind = 1 corresp. z = -1; ind = 2 to z = 1

    if (mean_minus == 0 && mean_plus == 0)
        kappa_lead = pi/2*(1 - 1/k^2)*cyc_dir*(amp_minus*(z_Kutta - 1/k) - amp_plus*(z_Kutta + 1/k));
        z_IC(ind) = z_Kutta + dt*sign(kappa_lead)*abs(kappa_lead)^(1/2)/(4*sqrt(2))*10^(1/4)*(sqrt(5/3) + sign(kappa_lead)*1i); 
    else
        kappa_lead = pi/2*(1 - 1/k^2)*(mean_minus*(z_Kutta - 1/k) - mean_plus*(z_Kutta + 1/k));
        z_IC(ind) = z_Kutta + sqrt(dt)*sign(kappa_lead)*abs(kappa_lead)^(1/2)/4*10^(1/4)*(sqrt(5/3) + sign(kappa_lead)*1i);
    end;

    zeta_IC(ind) = TJunctionPotential_Zeta_Func(z_IC(ind),k);
end;

disp(['vortex z-position is z = ' mat2str(z_IC) '; zeta = ' mat2str(zeta_IC)])
pause

[T,Z] = ode45(@(T,Z) RHSFunc(T,Z,k,mean_minus,amp_minus,mean_plus,amp_plus,tpe),tspan,z_IC,options);
disp(['size of Z is ']); size(Z)
disp('ode45 integration complete')

% calculate circulation values %
gamma_minus_vec = NaN*ones(1,length(T));
gamma_plus_vec = NaN*ones(1,length(T));
for tind = 1:length(T)
    u_minus = u_func(T(tind),mean_minus,amp_minus);
    u_plus = u_func(T(tind),mean_plus,amp_plus);
    
    z_minus = Z(tind,1);
    z_plus = Z(tind,2);
    
    if (tpe == 1)
        gamma_minus_vec(tind) = GammaCalcSingle(z_minus,u_minus,u_plus,k,-1);
        gamma_plus_vec(tind) = GammaCalcSingle(z_plus,u_minus,u_plus,k,1);
    else
        [gamma_minus,gamma_plus] = GammaCalc(z_minus,z_plus,u_minus,u_plus,k);
        gamma_minus_vec(tind) = gamma_minus;
        gamma_plus_vec(tind) = gamma_plus;
    end;
end;

% locate first time at which gamma_+ vanishes
loc = find(gamma_plus_vec(2:end).*gamma_plus_vec(1:end-1) < 0);
if (~isempty(loc))
    loc = loc(1);
    z_minus_loc = Z(loc,1);
    t_loc = T(loc);
    u_minus_loc = u_func(t_loc,mean_minus,amp_minus);
    u_plus_loc = u_func(t_loc,mean_plus,amp_plus);
    delta_loc = u_minus_loc - u_plus_loc;
    sigma_loc = u_minus_loc + u_plus_loc;
    gamma_minus_val1 = (delta_loc + sigma_loc/k)/(1 - 1/k^2)*abs(-1 - z_minus_loc)^2/imag(z_minus_loc);
    gamma_minus_val2 = (-delta_loc + sigma_loc/k)/(1 - 1/k^2)*abs(1 - z_minus_loc)^2/imag(z_minus_loc);
    disp(['gamma_- = ' num2str(gamma_minus_vec(loc)) ' should equal ' num2str(gamma_minus_val1) ' and ' num2str(gamma_minus_val2) '; gamma_+ = ' num2str(gamma_plus_vec(loc:loc+1))])
    zeta_minus_loc = TJunctionPotential_Zeta_Func(z_minus_loc,k);
    disp(['z_- = ' num2str(z_minus_loc) '; zeta_- = ' num2str(zeta_minus_loc)])
else
    disp(['gamma_+ is always nonzero'])
    loc = length(T);
end;

%%% plot solution in z-plane %%%

figure(cf)
subplot(5,3,[7 10])
hold on

plot(real(Z(:,1)),imag(Z(:,1)),'.-b','MarkerSize',16)
plot(real(Z(:,2)),imag(Z(:,2)),'.-r','MarkerSize',16)
hold off
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')

subplot(5,3,8)
plot(T,real(Z(:,1)),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($z_-$)','interpreter','latex')

subplot(5,3,9)
plot(T,imag(Z(:,1)),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($z_-$)','interpreter','latex')

subplot(5,3,11)
plot(T,real(Z(:,2)),'.-r','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($z_+$)','interpreter','latex')

subplot(5,3,12)
plot(T,imag(Z(:,2)),'.-r','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($z_+$)','interpreter','latex')

subplot(5,3,14)
plot(T,gamma_minus_vec,'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('$\gamma_-$','interpreter','latex')

subplot(5,3,15)
plot(T,gamma_plus_vec,'.-r','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('$\gamma_+$','interpreter','latex')

pause

% construct solution in zeta-plane
T = T(1:loc);
Z = Z(1:loc,:); % CHANGE LATER: TRUNCATE SOLUTION AT GAMMA_+ = 0 TIME POINT

zeta_mat = NaN*ones(size(Z));
for vind = 1:2 % loop over vortices
    for tind = 1:length(T) % loop over time
        zeta_mat(tind,vind) = TJunctionPotential_Zeta_Func(Z(tind,vind),k);
    end;
end;

%%% plot solution in zeta-plane %%%
%[gamma_minus,gamma_plus] = GammaCalc(z_minus,z_plus,u_minus,u_plus,k)

subplot(5,3,[1 4])
hold on

% construct lines showing T-junction boundary
line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
line([.5 2],[1 1],'LineWidth',3,'Color','k')
line([-2 2],[0 0],'LineWidth',3,'Color','k')

plot(real(zeta_mat(:,1)),imag(zeta_mat(:,1)),'.-b','MarkerSize',16)
plot(real(zeta_mat(:,2)),imag(zeta_mat(:,2)),'.-r','MarkerSize',16)

hold off

xlabel('Re($\zeta$)','interpreter','latex')
ylabel('Im($\zeta$)','interpreter','latex')

subplot(5,3,2)
plot(T,real(zeta_mat(:,1)),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($\zeta_-$)','interpreter','latex')

subplot(5,3,3)
plot(T,imag(zeta_mat(:,1)),'.-b','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($\zeta_-$)','interpreter','latex')

subplot(5,3,5)
plot(T,real(zeta_mat(:,2)),'.-r','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Re($\zeta_+$)','interpreter','latex')

subplot(5,3,6)
plot(T,imag(zeta_mat(:,2)),'.-r','MarkerSize',16)
xlabel('$t$','interpreter','latex')
ylabel('Im($\zeta_+$)','interpreter','latex')

function out = RHSFunc(t,z,k,mean_minus,amp_minus,mean_plus,amp_plus,tpe)
% right-hand side of vortex motion equation: dz/dt = ...

% important to initialize as a column vector (as is z)
out = NaN*ones(2,1); 

z_minus = z(1);
z_plus = z(2);
u_minus = u_func(t,mean_minus,amp_minus);
u_plus = u_func(t,mean_plus,amp_plus);

if (tpe == 1)
    % gamma_minus = 0; % CHANGE BACK LATER
    % gamma_minus = abs(-1 - z_minus)^2/imag(z_minus)*(-u_minus/(-1 + 1/k) + u_plus/(-1 - 1/k));
    % gamma_plus = 0;
    % gamma_plus = abs(1 - z_plus)^2/imag(z_plus)*(-u_minus/(1 + 1/k) + u_plus/(1 - 1/k));

    gamma_plus = 0;
    gamma_minus = GammaCalcSingle(z_minus,u_minus,u_plus,k,-1);
    rhs_minus = 1/abs(Fprime(z_minus,k))^2*(1/pi*(u_minus/(z_minus + 1/k) - u_plus/(z_minus - 1/k)) - 1i*gamma_plus/(2*pi)*(1/(z_minus - z_plus) - 1/(z_minus - conj(z_plus))) + 1i*gamma_minus/(2*pi*(z_minus - conj(z_minus))) + 1i*gamma_minus/(4*pi)*F2prime_over_Fprime(z_minus,k));

    gamma_minus = 0;
    gamma_plus = GammaCalcSingle(z_plus,u_minus,u_plus,k,1);
    rhs_plus = 1/abs(Fprime(z_plus,k))^2*(1/pi*(u_minus/(z_plus + 1/k) - u_plus/(z_plus - 1/k)) - 1i*gamma_minus/(2*pi)*(1/(z_plus - z_minus) - 1/(z_plus - conj(z_minus))) + 1i*gamma_plus/(2*pi*(z_plus - conj(z_plus))) + 1i*gamma_plus/(4*pi)*F2prime_over_Fprime(z_plus,k));
else
    [gamma_minus,gamma_plus] = GammaCalc(z_minus,z_plus,u_minus,u_plus,k);
    rhs_minus = 1/abs(Fprime(z_minus,k))^2*(1/pi*(u_minus/(z_minus + 1/k) - u_plus/(z_minus - 1/k)) - 1i*gamma_plus/(2*pi)*(1/(z_minus - z_plus) - 1/(z_minus - conj(z_plus))) + 1i*gamma_minus/(2*pi*(z_minus - conj(z_minus))) + 1i*gamma_minus/(4*pi)*F2prime_over_Fprime(z_minus,k));
    rhs_plus = 1/abs(Fprime(z_plus,k))^2*(1/pi*(u_minus/(z_plus + 1/k) - u_plus/(z_plus - 1/k)) - 1i*gamma_minus/(2*pi)*(1/(z_plus - z_minus) - 1/(z_plus - conj(z_minus))) + 1i*gamma_plus/(2*pi*(z_plus - conj(z_plus))) + 1i*gamma_plus/(4*pi)*F2prime_over_Fprime(z_plus,k));
end;

out(1) = conj(rhs_minus);
out(2) = conj(rhs_plus);

function [gamma_minus,gamma_plus] = GammaCalc(z_minus,z_plus,u_minus,u_plus,k)

vel_delta = u_minus - u_plus;
vel_sigma = u_minus + u_plus;

% solve system of 2 equations and two unknowns for vortex strengths
%Mat(1,1) = imag(z_minus)/abs(-1-z_minus)^2;
%Mat(1,2) = imag(z_plus)/abs(-1-z_plus)^2;
%Mat(2,1) = imag(z_minus)/abs(1-z_minus)^2;
%Mat(2,2) = imag(z_plus)/abs(1-z_plus)^2;
%rhs_vec = 1/(1 - 1/k^2)*[vel_delta + vel_sigma/k; -vel_delta + vel_sigma/k];
%gamma_vec = inv(Mat)*rhs_vec;
%gamma_minus = gamma_vec(1);
%gamma_plus = gamma_vec(2);

% use explicit solution for above system of equations
denom = abs(1 + z_plus)^2*abs(1 - z_minus)^2 - abs(1 + z_minus)^2*abs(1 - z_plus)^2;
gamma_minus = abs(1 - z_minus^2)^2/(imag(z_minus)*denom*(1 - 1/k^2))*(vel_delta*(abs(1 + z_plus)^2 + abs(1 - z_plus)^2) + vel_sigma/k*(abs(1 + z_plus)^2 - abs(1 - z_plus)^2));
gamma_plus = -abs(1 - z_plus^2)^2/(imag(z_plus)*denom*(1 - 1/k^2))*(vel_delta*(abs(1 + z_minus)^2 + abs(1 - z_minus)^2) + vel_sigma/k*(abs(1 + z_minus)^2 - abs(1 - z_minus)^2));

function gamma_out = GammaCalcSingle(z,u_minus,u_plus,k,z_Kutta)

gamma_out = abs(z_Kutta - z)^2/imag(z)*(-u_minus/(z_Kutta + 1/k) + u_plus/(z_Kutta - 1/k));

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