function classpractice51()
  
clc
clear all

disp('----------------------------------------------------------------------');
disp('-------------------------class practice 5.1----------------------------');
disp('-------------------by Amber-GaoQi on 23/8/26---------------------');
disp('-----------------------------------------------------------------------');
disp('0-SWF');
disp('1-tx-ITL with fixed S');
disp('2-self-optimizing tx-ITL');
disp('3-FLG');
disp('4-FFR+IFR3');
disp('5-FFR+SWF');




while 1
    mOPTION=input('choose 0,1,2,3,6,7: ');
    if mOPTION==0||mOPTION==1||mOPTION==2||mOPTION==3
        break
    end
end

OPTIONS(1)=mOPTION;
if mOPTION==1||mOPTION==2
    S=input('the transmit interference temperature limit: ');
    OPTIONS(2)=S;
    if mOPTION==2
        beta=input('the beta: ');
        OPTIONS(3)=beta;
    end
end

B=input('number of cells or BSs: ');
Um=input('the average number of users in a cell');
Nf=input('number of frequency cells');
Nsim=input('number of random simulation');
aveSNRdB=input('average SNR per user in dB');

aveSNR=10.^(aveSNRdB/10);
P=1;
sigma=P./sqrt(aveSNR);
fprintf('\n');

%Nsim random simulation for average capacity
for n=1:Nsim
    %for each ransom simulation
    %generate the number of users in each cell
    for b=1:B
        U(b)=poissrnd(Um);
    end
    
    %generate the channels
    for b_=1:B
        for b =1:B
            for u=1:U(b)
                H(b_,b,u,1:Nf)=(randn(1,Nf)+1i*randn(1,Nf))/sqrt(2);
            end
        end
    end
    
    %power allocation with different method
    popt=sWF1(H,sigma,P,Nf,B,U,OPTIONS); 
    
    %compute the capacity
    for b=1:B
        for u=1:U(b)
            C_user(b,u)=0;
            
            %calculate the interference
            for f=1:Nf
                I(f)=0;
                
                for b_=1:b
                    for u_=1:U(b_)
                        if ~((b_==b)&&(u_==u))
                            I(f)=I(f)+popt(b_,u_,f)*abs(H(b_,b,u,f))^2;
                        end
                    end
                end
                
                C_user(b,u)=C_user(b,u)+log2(1+popt(b,u,f)*abs(H(b,b,u,f))^2/(I(f)+sigma^2));
            end
        end
        
        if U(b)>=1
            C_cell(b)=sum(C_user(b,1:U(b)))
         else
             C_cell(b)=0
        end
        
    end
    
    %the capacity of the network
    C_net(n)=sum(C_cell);
    
    if mod(n,100)==0
        fprintf('\n');
        fprintf('%d simulation sompleted\n', n);
    end
    

%the average network capacity


fprintf('\n');
fprintf('the average network capacity is %g', ave_C_Net);

        
    
end