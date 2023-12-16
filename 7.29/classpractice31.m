function classpractice31()

clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 3.1(19 cells )----------------------------')
disp('-------------------by Amber-GaoQi on 23/8/12---------------------')
disp('-----------------------------------------------------------------------')
disp('1-IfR1')
disp('3-IfR3')

which_method=input('choose the number: ');
fprintf('\n')


%Nsim=input('the number of random simulations: ');

switch which_method
    case 0 %IFR
        reuse_factor=input('reuse factor, 1 or 3: ');
        Nc=input('the number of frequency channels: ');
        Um=input('average number of users in each cell: ');
        snrdB=input('average SNR for each user in dB: ');
        snr=10^(snrdB/10);
        P=1;
        sigma=sqrt(P/snr);
        
        if reuse_factor==1
            
            for c=1:19
                f(c,1:Nc)=1:Nc; %frequency matrix 
            end
            r=0.7;                  %rho matrix
            
            rho = [
                1 r r r r r r r r r r r r r r r r r;...
                r 1 r r r r r r r r r r r r r r r r;...
                r r 1 r r r r r r r r r r r r r r r;...
                r r r 1 r r r r r r r r r r r r r r;...
                r r r r 1 r r r r r r r r r r r r r;...
                r r r r r 1 r r r r r r r r r r r r;...
                r r r r r r 1 r r r r r r r r r r r;...
                r r r r r r r 1 r r r r r r r r r r;...
                r r r r r r r r 1 r r r r r r r r r;...
                r r r r r r r r r 1 r r r r r r r r;...
                r r r r r r r r r r 1 r r r r r r r;...
                r r r r r r r r r r r 1 r r r r r r;...
                r r r r r r r r r r r r 1 r r r r r;...
                r r r r r r r r r r r r r 1 r r r r;...
                r r r r r r r r r r r r r r 1 r r r;...
                r r r r r r r r r r r r r r r 1 r r;...
                r r r r r r r r r r r r r r r r 1 r;...
                r r r r r r r r r r r r r r r r r 1];

        else
            Nc1=round(Nc/3);
            Nc2=Nc1;
            Nc3 = Nc-Nc1-Nc2;
            f(1,1:Nc1)=1:Nc1;
            f(2,1:Nc2)=Nc1+1:Nc1+Nc2;
            f(3,1:Nc3)=Nc1+Nc2+1:Nc;
            f(4, : )=f(2, : );
            f(5, : )=f(3, : );
            f(6, : )=f(2, : );
            f(7, : )=f(3, : );
            f(8, : )=f(1, : );
            f(9, : )=f(3, : );
            f(10, : )=f(1, : );
            f(11, : )=f(2, : );
            f(12, : )=f(1, : );
            f(13, : )=f(3, : );
            f(14, : )=f(1, : );
            f(15, : )=f(2, : );
            f(16, : )=f(1, : );
            f(17, : )=f(3, : );
            f(18, : )=f(1, : );
            f(19, : )=f(2, : );
            r=0.7;
          %  rho=[1 0 0 0 0 0 0 r 0 r 0 r 0 r 0 r 0 r 0 ;...
             %        0 1 0 r 0 r 0 0 0 0 r 0 0 0 r 0 0 0 r;...
                %     0 0 1 0 r 0 r 0 r 0 0 0 r 0 0 0 r 0 0;...
                   %  0 r 0 1 0 r 0 0 0 0 r 0 0 0 r 0 0 0 r;...
                   %0 0 r 0 1 0 r 0 r 0 0 0 r 0 0 0 r 0 0;...
                    % 0 r 0 r 0 1 0 0 0 0 r 0 0 0 r 0 0 0 r;...
                     %0 0 r 0 r 0 1 0 r 0 0 0 r 0 0 0 r 0 0];
                 
            rho = [
                1 0 0 0 0 0 0 r 0 r 0 r 0 r 0 r 0 r 0;
                0 1 0 r 0 r 0 0 0 0 r 0 0 0 r 0 0 0 r;
                0 0 1 0 r 0 r 0 r 0 0 0 r 0 0 0 r 0 0;
                0 r 0 1 0 r 0 0 0 0 r 0 0 0 r 0 0 0 r;
                0 0 r 0 1 0 r 0 r 0 0 0 r 0 0 0 r 0 0;
                0 r 0 r 0 1 0 0 0 0 r 0 0 0 r 0 0 0 r;
                0 0 r 0 r 0 1 0 r 0 0 0 r 0 0 0 r 0 0;
                r 0 0 0 0 0 0 1 0 r 0 r 0 r 0 r 0 r 0;
                0 0 r 0 r 0 r 0 r 0 0 0 0 r 0 0 0 r 0;
                r r 0 r 0 r 0 0 0 0 r 0 0 0 r 0 0 0 r;
                0 0 0 0 0 0 0 r 0 1 0 r 0 r 0 r 0 r 0;
                r 0 r r 0 r 0 0 0 r 0 0 0 0 r 0 0 0 r;
                0 0 0 0 r 0 r 0 r 0 r 0 1 0 r 0 0 0 r;
                r r 0 r 0 r 0 0 0 r 0 0 r 0 r 0 0 0 r;
                0 0 0 0 0 0 0 r 0 r 0 r 0 0 1 0 r 0 r;
                r 0 r r 0 r 0 0 0 r 0 0 0 0 r 0 0 0 r;
                0 0 0 0 r 0 r 0 r 0 r 0 r 0 r 0 1 0 r;
                r r 0 r 0 r 0 0 0 r 0 0 r 0 r 0 0 0 r;
                0 0 0 0 0 0 0 r 0 r 0 r 0 r 0 r 0 r 1];

        end

    %setting the number of active users for all the cells
    for c=1:19
        U(c)=poissrnd(Um);
        Nchannel(c)=length( find( f(c , : )>0));%check
        if U(c)>Nchannel(c)
            fprintf('we cannot server all the %d users in cell %d\n',Nchannel(c),c)
        end
    end
    
    %generate fading channels for each serving user from all base stations
    for c1=1:19
        for c2=1:19
            for u=1:min(Nchannel(c2),U(c2))
                h(c1,c2,u)=(randn+1i*randn)/sqrt(2); %from cell c1 to user u in cell c2
            end
        end
    end
    
    %compute the user capacity for each cell
    for c=1:19
        for u=1:min(Nchannel(c),U(c))
            I(c,u)=0;
            for c_=1:19
                if rho(c_,c)~=0&&c_~=c
                    if U(c_)>=u&&Nchannel(c_)>=u
                         I(c,u)=I(c,u)+rho(c_,c)*P*abs(h(c_,c,u))^2;
                    end
                end
            end
            sinr(c,u)=P*abs(h(c,c,u))^2/(I(c,u)+sigma^2);
            capacity(c,u)=log2(1+sinr(c,u));
        end
        cell_capacity(c)=sum(capacity(c, : ));
    end
    total_net=sum(cell_capacity);
    fprintf('capacity is %g\n',total_net);
   
  
  case 1
      
  case 2
          
   
end
end

    