function fval=fun_waterlevel(x,sigma,H,P)

Nf=length(H);
fval=0;
for f=1:Nf
    fval=fval+max([0 x-sigma^2/abs(H(f))^2]);
end

fval=fval-P;

end