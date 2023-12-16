function popt=sWF1(H,sigma,P,Nf,B,U,OPTIONS)

if OPTIONS(1)==0
    mOPTION=0;
elseif OPTIONS(1)==1
    mOPTION=1;
    S=OPTIONS(2);
elseif OPTIONS(1)==2
    mOPTION=2;
    S=OPTIONS(2);
    beta=OPTIONS(3);
elseif OPTIONS(1)==3
    mOPTIONS=3;
end


Maxiter=200;
delta=le-4;

for b=1:B
    for u=1:U(b)
        ptmp=rand(1,Nf);
        ptmp=P*ptmp/sum(ptmp);
        p(b,u,1:Nf)=ptmp;
        clear ptmp;
    end
end

n_iter=1;
while 1
    p_old=p;
    
    for b=1:B
        for u=1:U(b)
            %focusing a certain user
            if mOPTION==0%SWF
                ptmp=findp_wI(p,H,sigma,P,Nf,B,U,b,u);
            elseif mOPTION==1%txITL
                ptmp=finp_txITL(p,H,sigma,P,Nf,B,U,b,u,S);
            elseif mOPTION==2%self-txITL
                ptmp=findp_txITL_beta(p,H,sigma,P,Nf,B,U,b,u,S,beta);
            elseif mOPTION==3%FLG
                ptmp=findp_FLG(p,H,sigma,P,Nf,B,U,b,u)
            end
            p(b,u,1:Nf)=ptmp;
            clear ptmp;
        end
    end
    
    p_new=p;
    n_iter=n_iter+1;
    
    if n_iter>Maxiter||sum(sum(sum((p_new-p_old).^2)))<delta
        break
    end
    
end

end