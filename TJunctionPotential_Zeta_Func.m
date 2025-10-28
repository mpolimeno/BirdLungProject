function out = TJunctionPotential_Zeta_Func(z,k)

%if (imag(z) ~= 0 || abs(real(z)) < 1/k)
if (imag(z) ~= 0)
    s = @(t) t*z;
    ds = @(t) z;
else
    if (real(z) >= 0)
        s = @(t) z/2 * (1+exp(1i*pi*(1-t)));
        ds = @(t) -1i*pi*z/2 * exp(1i*pi*(1-t));
    else
        s = @(t) abs(z)/2 * (-1 + exp(1i*pi*t)); % MATTEO: exp(1i*pi*t)=(-1)^t
        ds = @(t) 1i*pi*abs(z)/2 * exp(1i*pi*t);
    end;
end;
%out = 5/pi*quadgk(@(t) sqrt(1-s(t)).*sqrt(1+s(t))./(1-(k*s(t)).^2).*ds(t),0,1);
out = 5/pi*quadgk(@(t) sqrt(1-s(t).^2)./(1-(k*s(t)).^2).*ds(t),0,1); %%Matteo: this uses the Gauss-Kronrod quadrature to evaluate that integral 
% (which is the conformal map from the complex half plane to the T junction)