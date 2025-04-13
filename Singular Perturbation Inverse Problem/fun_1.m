function val=fun_1(x,t,epsilon)
val=exp(-t).*(exp(-1./epsilon)+(1-exp(-1./epsilon)).*x-exp(-(1-x)./epsilon));
end