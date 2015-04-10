%******************************************* 
%   LATTICE BOLTZMANN MODEL 2D             *
%******************************************* 
% numerotation des vitesses 
%     6  2  5 
%      \ | / 
%       \|/ 
%   3----9----1 
%       /|\ 
%      / | \ 
%     7  4  8
%*******************************************
%              Collision D2Q9 cs2 reglable *
%                                          * 
%*******************************************
%   nx, ny : TAILLE DU RESEAU
%   N : DISTRIBUTIONS SUR LE RESEAU
%   PHY : COEFFICIENT PRE_CALCULES
%   
%   PHY(1,1)=\rho/3
%   PHY(1,2)=(4-6 cs2) \rho 
%   PHY(1,3)=(4-9 cs2) \rho
%   PHY(1,4)=3/ \rho
%   PHY(1,5)=1/ \rho
%   PHY(1,6)=LAMBDA_S/9
%   PHY(1,7)=LAMBDA_S/18
%   PHY(2,1)=LAMBDA_E/36
%   PHY(2,2)=LAMBDA_S/36
%   PHY(2,3)=LAMBDA_X/12
%   PHY(2,4)=LAMBDA_NU/4
%   PHY(2,5)=LAMBDA_E/18
%   PHY(2,6)=LAMBDA_E/9
%   PHY(2,7)=LAMBDA_X/6
%  DESCRIP-END.

%  9th FEB 2015: all control=OK
% constant density, no mass loss
% relative error stable at 1e-3


close('all');clear('all');

% init Lattice main values
nx=5; ny=31;  % lattice size
CS2=1.0/3.0; % speed of sound
DENSITE = 1.0;
LAMBDA = 0.2;
SIGMA = LAMBDA;
LAMBDAE = 0.2;
LADEUX = 0.2;
TAU=(4.0-2.0*LAMBDA)/(2.0+LAMBDA*(3.0*LADEUX-1.0));
t1=(3.0-5.0*CS2)/3.0; 
t2=CS2/3.0; 
t3=CS2/12.0; 

PHY(1,1)=DENSITE/3.0;
PHY(1,2)=(4.0-6.0*CS2)*DENSITE;
PHY(1,3)=(4.0-9.0*CS2)*DENSITE;
PHY(1,4)=3.0/DENSITE;
PHY(1,5)=1.0/DENSITE;
PHY(1,6)=SIGMA/9.0;
PHY(1,7)=SIGMA/18.0;
PHY(2,1)=LAMBDAE/36.0;
PHY(2,2)=SIGMA/36.0;
PHY(2,3)=TAU/12.0;
PHY(2,4)=LAMBDA/4.0;
PHY(2,5)=LAMBDAE/18.0;
PHY(2,6)=LAMBDAE/9.0;
PHY(2,7)=TAU/6.0;

TN = [ t2,t2,t2,t2,  t3,t3,t3,t3, t1];%weight of different speeds
CX = [    1,  0, -1,  0,    1,  -1,  -1,   1,0];
CY = [    0,  1,  0, -1,    1,   1,  -1,  -1,0];
opp = [ 1, 4, 5, 2,  3, 8, 9,  6,  7]; 
BOUND=zeros([nx ny]);
max_t = 4000;

% init Population
N=zeros([nx ny 9]);
F=zeros([nx ny 9]);
for IX = 1: nx
    BOUND(IX,1)=1; % Poiseuil flow
    BOUND(IX,ny)=1; % Poiseuil flow
    
         for IY = 1: ny
            for dir = 1: 9
              N(IX,IY,dir)=TN(dir)*DENSITE;  
            end
         end
end

 for j=1:ny 
        for i=1:nx
            %add space matricies for plotting
            x(j,i)=i;
            y(j,i)=j;
        end
    end


       

avu=1; 
prevavu=1; ts=0; deltaUx=1e-5; deltaUy=0; % 1e-10
DENSITY=sum(N,3);
UX=zeros([nx ny]);
UY=zeros([nx ny]);
numactivenodes=sum(sum(1-BOUND));
while (ts<max_t) & 1e-7<abs((prevavu-avu)/avu) | ts<100

    % collision 
    for IX = 1: nx
         for IY = 1: ny
            
            B0=N(IX,IY,9);
            B1=N(IX,IY,1);
            B2=N(IX,IY,2);
            B3=N(IX,IY,3);
            B4=N(IX,IY,4);
            B5=N(IX,IY,5);
            B6=N(IX,IY,6);
            B7=N(IX,IY,7);
            B8=N(IX,IY,8);
             if ( BOUND(IX,IY)==0) % not a boundary
            RO=B0+(B1+B2)+(B3+B4)+(B5+B6)+(B7+B8);
            JX=B1-B3+B5-B6-B7+B8;
            JY=B2-B4+B5+B6-B7-B8;
            UX(IX,IY)=JX;
            UY(IX,IY)=JY;
            DENSITY(IX,IY)=RO;
	    E=4.0*B0 +(B1+B2+B3+B4)-2.0*(B5+B6+B7+B8)-PHY(1,2)*RO...
            +PHY(1,4)*(JX*JX+JY*JY);
	    H=4.0*B0 -2.0*(B1+B2+B3+B4)+(B5+B6+B7+B8)-PHY(1,3)*RO...
            +PHY(1,4)*(JX*JX+JY*JY);
	    XX=2.0*(B1-B3)-B5+B6+B7-B8-JX;
	    XY=2.0*(B2-B4)-(B5+B6)+(B7+B8)-JY;
	    SX=B1-B2+B3-B4-PHY(1,5)*(JX*JX-JY*JY);
	    SY=B5-B6+B7-B8-PHY(1,5)*(JX * JY);

	    B0=B0-PHY(2,6)*E-PHY(1,6)*H;
	    B1=B1-PHY(2,1)*E+PHY(1,7)*H-PHY(2,7)*XX-PHY(2,4)*SX;
	    B2=B2-PHY(2,1)*E-PHY(2,7)*XY+PHY(2,4)*SX+PHY(1,7)*H;
	    B3=B3-PHY(2,1)*E+PHY(2,7)*XX-PHY(2,4)*SX+PHY(1,7)*H;
	    B4=B4-PHY(2,1)*E+PHY(2,7)*XY+PHY(2,4)*SX+PHY(1,7)*H;
	    B5=B5+PHY(2,5)*E+PHY(2,3)*(XX+XY)-PHY(2,4)*SY-PHY(2,2)*H;
	    B6=B6+PHY(2,5)*E+PHY(2,3)*(XY-XX)+PHY(2,4)*SY-PHY(2,2)*H;
	    B7=B7+PHY(2,5)*E-PHY(2,3)*(XX+XY)-PHY(2,4)*SY-PHY(2,2)*H;	    	    
	    B8=B8+PHY(2,5)*E-PHY(2,3)*(XY-XX)+PHY(2,4)*SY-PHY(2,2)*H;
            F(IX,IY,9)=B0;
            F(IX,IY,1)=B1;
            F(IX,IY,2)=B2;
            F(IX,IY,3)=B3;
            F(IX,IY,4)=B4;
            F(IX,IY,5)=B5;
            F(IX,IY,6)=B6;
            F(IX,IY,7)=B7;
            F(IX,IY,8)=B8;
             else % bounce back on boundaries -> revert populations
            F(IX,IY,9)=B0;
            F(IX,IY,1)=B3;
            F(IX,IY,2)=B4;
            F(IX,IY,3)=B1;
            F(IX,IY,4)=B2;
            F(IX,IY,5)=B7;
            F(IX,IY,6)=B8;
            F(IX,IY,7)=B5;
            F(IX,IY,8)=B6; 
            %UX(IX,IY)=0.0;
            %UY(IX,IY)=0.0;
            DENSITY(IX,IY)=DENSITE;
             end
         end
    end
     %Forcing 
    %
    for IX = 1: nx
         for IY = 1: ny
             if ( BOUND(IX,IY)==0) % not a boundary
                for dir=1:9
                    F(IX,IY,dir)=F(IX,IY,dir)+ TN(dir)* DENSITE*(deltaUx*...
                        CX(dir)+deltaUy*CY(dir))/CS2;
                end
             end
         end
    end
    
  
% STREAMING STEP 
for i=1:9 
    N(:,:,i ) = circshift(F(:,:,i ), [CX(i),CY(i),0]); 
end


    prevavu=avu;
    avu=sum(sum(UX))/numactivenodes; % test convergence
    ts=ts+1;
end

figure;colormap(gray(2));image(2-BOUND');hold on;
quiver(2:nx,1:ny,UX(2:nx,:)',UY(2:nx,:)');
title(['Flow field after ',num2str(ts),'\deltat']);xlabel('x');ylabel('y');

width=ny-2;

figure
%Model results as red circles
%plot(y(:,1)-1.5-width/2,UX(:,1),'ro');
plot(y([2:ny-1],1)-1.5-width/2,UX(1,[2:ny-1])+deltaUx/2,'ro');
hold on

%Poiseuille velocity profile as blue line

nu=1/3*(1.0/LAMBDA-1/2);
g=deltaUx;
plot(y(:,1)-1.5-width/2,g/(2*nu)*((width/2)^2-(y(:,1)-1.5-width/2).^2)) ;
title(['Velocity Profil after ',num2str(ts),'\deltat']);xlabel('x');ylabel('y');

%print error
error = zeros(ny);
for i= 2:ny-1
   error(i)= (g/(2*nu)*( (width/2)^2-(y(i,1)-1.5-width/2).^2)-(UX(1,i)...
       +deltaUx/2))/(UX(1,i)+deltaUx/2);
end
figure
plot(y(2:ny-1,1)-1.5-width/2,error(2:ny-1),'ro') ;
title(['Relative Error in Velocity after ',num2str(ts),'\deltat']);xlabel('x');ylabel('y');

