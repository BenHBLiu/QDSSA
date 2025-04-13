function val=fun_ue(x,t,epsilon)
val=(1-exp(-t)).*((exp(-x./sqrt(epsilon))+exp(-(1-x)./sqrt(epsilon)))./(1+exp(-1./sqrt(epsilon)))-cos(pi.*x).^2);
end