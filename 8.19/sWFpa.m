function p=sWFpa(H,sigma,P)
%compute the power allocation

[Nf,U,~]=size(H);
Niter_max=1e2;
d_threshold=1e-4;
Nrun=1;
p=rand(Nf,U)*P/Nf;
p_old=p;

while 1
    for u=1:U
        for f=1:Nf
            I(f)=sigma^2;
            for u_=1:U
                if u_~=u
                    I(f)=I(f)+p(f,u_)*abs(H(f,u_,u))^2;
                end
            end
            c(f)=I(f)/abs(H(f,u,u))^2;
        end
        ptmp=wfpa_c(c,P);
        p(:,u)=ptmp';
        clear ptmp;
    end
    Nrun=Nrun+1;
    p_new=p;
    if sum(sum(abs(p_new-p_old).^2))<d_threshold | Nrun>Niter_max
        return;
    end
end
end