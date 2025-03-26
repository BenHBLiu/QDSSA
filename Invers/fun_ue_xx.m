function val=fun_ue_xx(x,t,epsilon)
val=-(exp(-t)-1).*((exp(-x./epsilon.^(1./2))./epsilon + exp((x-1)./epsilon.^(1./2))./epsilon)./(exp(-1./epsilon.^(1./2))+1)-2.*pi.^2.*sin(pi.*x).^2+2.*pi.^2.*cos(pi.*x).^2);
end