function [X1,Y1,Z1,X2,Y2,Z2] = TJunctionPotential_3

options = optimset('display','off');

k = sqrt(5);

d_space = .01;

% horizontal branch of pipe first; so iterate using vertical lines
xi_v = -2:d_space:2; % variable along horizontal branch (x)
eta_v = d_space:d_space:(1-d_space); % variable along pipe diameter (y)
[X1,Y1] = meshgrid(xi_v,eta_v);
Z1 = NaN*ones(size(X1));
for xind = 1:length(xi_v)
    disp(['horizontal branch: xi = ' num2str(xi_v(xind))]);

    zg = [0 1]; % at start of each line, restart guess
    
    for ind = 1:length(eta_v)
        zeta = xi_v(xind) + 1i*eta_v(ind);
        [z,~,eflag] = fsolve(@(Z) RootFunc(Z,zeta,k),zg,options);
        if (eflag > 0)
            %z_vec(ind) = z(1) + 1i*z(2);
            zg = z;
            Z1(ind,xind) = z(1) + 1i*z(2);
        end;
        %disp(['eta = ' num2str(eta_v(ind)) '; eflag = ' num2str(eflag) '; z = ' num2str(z_vec(ind))]);
    end;
end;

% vertical branch of pipe for Im(zeta) > 1; iterate using horizontal lines
eta_v = (-.5+d_space):d_space:(.5-d_space); % variable along pipe diameter (x)
xi_v = (1+d_space):d_space:3; % variable along vertical branch (y)
[X2,Y2] = meshgrid(eta_v,xi_v);
Z2 = NaN*ones(size(X2));
for xind = 1:length(xi_v)
    disp(['vertical branch: xi = ' num2str(xi_v(xind))]);
    zg = [0 1]; % at start of each line, restart guess
    for ind = 1:length(eta_v)
        zeta = eta_v(ind) + 1i*xi_v(xind);
        [z,~,eflag] = fsolve(@(Z) RootFunc(Z,zeta,k),zg,options);
        if (eflag > 0)
            %z_vec(ind) = z(1) + 1i*z(2);
            zg = z;
            Z2(xind,ind) = z(1) + 1i*z(2);
        end;
    end;
end;

function out = RootFunc(zvec,zeta_val,k)

z = zvec(1) + 1i*zvec(2);
zeta = TJunctionPotential_Zeta_Func(z,k);
out(1) = real(zeta - zeta_val);
out(2) = imag(zeta - zeta_val);