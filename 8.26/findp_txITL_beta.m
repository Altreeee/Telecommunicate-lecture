function popt=findp_txITL_beta(p,H,sigma,P,Nf,B,U,b,u,S,beta)


%interference
for f=1:Nf
    I(f)=0;
    for b_=1:N
        for u_=1:U(b_)
            if ~(b_==b0)&&(u_==u0)
                I(f)=I(f)+p(b_,u_,f)*abs(H(b_,b0,u0,f))^2;
            end
        end
    end
end

Htmp(1:Nf)=H(b0,b0,u0,1:Nf);
wl=fzero(@(w) fPcon_txITL_beta(w,Htmp,sigma,I,P,S,beta,p,u0,b0),P/Nf);

%average  interference
Ave_I=0;
N_p=0;
for f=1:Nf
    if p(b0,u0,f)>0
        Ave_I=Ave_I+I(f);
        N_p=N_p+1;
    end
end
Ave_I=Ave_I/N_p;
S_=max(0,Ave_I*beta)

for f=1:Nf
    c(f)=(I(f)+sigma^2)/abs(H(b0,b0,u0,f))^2;
    
    if c(f)<=S_
        popt(f)=max(0,wl-c(f));
    else
        popt(f)=0;
    end
    
end




end