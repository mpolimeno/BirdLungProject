function PipeModel_ODE_Vort_1

% need to handle initial time step ~ sqrt(t) (assuming mean is nonzero)
% need to loop over half-cycles to re-initialize vortices
% insert parameter values

k = sqrt(5);



function out = RHS_func(t,X,kL,Ls,ds,ks,d,z_Kutta,k)

out = NaN*ones(3,1);

u_s = X(1);
z_A = X(2);
z_B = X(3); 

ds_tilde = ds - 2*d;
Ls_tilde = Ls - 2*d;

u = u_func(t);
dot_u = dot_u_func(t);

[pA2A1,pA3A1] = PressureDiff(u,u-u_s,k,z_A,d);
[pB2B1,pB3B1] = PressureDiff(-u,-u_s,k,z_B,d);

pdiff_total = pA3A1 - pA2A1 + pB3B1 - pB2B1;

us_coeff = (kL*Ls_tilde + ks*ds_tilde)/(Ls_tilde + ds_tilde);
dot_u_coeff = ds_tilde/(Ls_tilde + ds_tilde);
u_coeff = ks*dot_u_coeff;
p_coeff = 1/(Ls_tilde + ds_tilde);

out(1) = -us_coeff*u_s + dot_u_coeff*dot_u + u_coeff*u + p_coeff*pdiff_total;
out(2) = VortexVelocityTJunc(z_A,u,u-u_s,k,z_Kutta);
out(3) = VortexVelocityTJunc(z_B,-u,-u_s,k,z_Kutta);

function out = u_func(t,A)

out = A*sin(t);

function out = dot_u_func(t,A)

out = A*cos(t);

function out = VortexVelocityTJunc(z,u_minus,u_plus,k,z_Kutta)

gamma = GammaCalcSingle(z,u_minus,u_plus,k,z_Kutta);
conj_vel = 1/abs(Fprime(z,k))^2*(1/pi*(u_minus/(z + 1/k) - u_plus/(z - 1/k)) + 1i*gamma/(2*pi*(z - conj(z))) + 1i*gamma/(4*pi)*F2prime_over_Fprime(z,k));
out = conj(conj_vel);

function out = ComplexVelocity(Z,k,u_minus,u_plus,gamma,z_vort) % velocity u-i*v = dw/dzeta

% term1 is "velocity" term (no vortex)
term1 = (-u_minus*(Z - 1/k) + u_plus*(Z + 1/k))./sqrt(1 - Z.^2);

% term2 incluces the vortex
term2 = -1i*gamma/(2*k^2)*(1 - (k*Z).^2)./sqrt(1 - Z.^2).*(1./(Z - z_vort) - 1./(Z - conj(z_vort)));

out = term1 + term2;

function [p_X2X1,p_X3X1] = PressureDiff(u_minus,u_plus,k,z_vort,d)
% think of "X" as "A" or "B" junction

load('TJunctionPotential_3-1a.mat');

gamma = GammaCalcSingle(z,u_minus,u_plus,k,z_Kutta);

CVel1 = ComplexVelocity(Z1,k,u_minus,u_plus,gamma,z_vort);
CVel2 = ComplexVelocity(Z2,k,u_minus,u_plus,gamma,z_vort);

X1_loc = find(ReZeta1(end,:) == (-0.5 - d));
X2_loc = find(ReZeta1(end,:) == (0.5 + d));
X3_loc = find(ImZeta2(:,end) == (1 + d));
SpeedSq1 = mean(abs(CVel1).^2,1);
SpeedSq2 = mean(abs(CVel2).^2,2);
p_X2X1 = 1/2*(SpeedSq1(X1_loc) - SpeedSq1(X2_loc));
p_X3X1 = 1/2*(SpeedSq1(X1_loc) - SpeedSq2(X3_loc));

function out = Fprime(z,k)
% derivative of conformal map from upper half-plane to T-junction

out = 5/pi*sqrt(1 - z.^2)./(1 - (k*z).^2);

function gamma_out = GammaCalcSingle(z,u_minus,u_plus,k,z_Kutta)

gamma_out = abs(z_Kutta - z)^2/imag(z)*(-u_minus/(z_Kutta + 1/k) + u_plus/(z_Kutta - 1/k));

function out = F2prime(z,k)
% second derivative of conformal map from upper half-plane to T-junction

out = -5/pi*z./sqrt(1 - z.^2).*(1 + k^2*(z.^2 - 2))./(1 - (k*z).^2).^2;

function out = F2prime_over_Fprime(z,k)

out = -z/(1 - z^2)*(1 + k^2*(z^2 - 2))/(1 - (k*z)^2);