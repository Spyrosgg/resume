---
layout: post
title:  "My Dissertation"
date:   2018-10-28 15:00:00 +0100
---

## Program Velocity Profile [VePro]
The program VePro, is using a  linearized fully parabolic version of the  depth-averaged  flow  equations, calculates an analytical solution for the transverse distribution of longitudinal velocity in  rectangular channels  and  as well as in channels of compound cross-section.
The program  also calculates the  discharge with the Analytical method and with the Manning's formula. Here  we give an  example  for a channel with two  cross-sections.  

### MATLAB code
    %_____________________________# INPUT #________________________________  
    % Insert the sum of the cross sections 
    cs=2;
    
    % Insert the heights hi=[h1,h2,h3....hv]
    hi=[0.025,0.05];
    
    % Insert the width of each cross section Dy=[Dy1,Dy2,Dy3....Dyv]
    Dy=[0.19,0.37];
    
    % Insert the manning number ni=[n1,n2,n3....nv]
    ni=[0.0116,0.0116];
    
    % Insert the discretization accuracy
    Sum=1000;
    
    % Insert the slope
    s=0.00047;
      
    % Insert the Gravity (m/s^2)
    g=9.81;
    
    % Insert the cv=[cv1,cv2,cv3....cvv]
    Cv=[0.20,0.20];
  
    >__________________________________## BEGIN ##_____________________________  

    %Setting the matrices where we will save the parameters
     h=[hi,0];
     n=[ni,0];
     a=zeros(1,cs);
     b=zeros(1,cs);
     Cf=zeros(1,cs);
     gt=zeros(1,cs);
     H=zeros(1,cs+1);
     Ai=zeros(1,cs);
     Bi=zeros(1,cs);
     Ci=zeros(1,cs);
     Di=zeros(1,cs);
     Gt1=zeros(1,cs);
     Gt2=zeros(1,cs);
     Ca=zeros(1,cs);
     Cb=zeros(1,cs);
      
     % Calculate the �, b, � parameters for each cross-section  
     
     for k=1:cs
         
         Cf(k)=n(k)^2*g/(h(k)^(1/3));    
         a(k)=sqrt(2*sqrt(Cf(k))/(Cv(k)*h(k)^2));    
         b(k)=-a(k)^2*g*(h(k))*s/Cf(k);
         gt(k)=-b(k)/a(k)^2;
 
         H(k+1)=(n(k+1)/n(k))*(h(k+1)/h(k))^(5/6);
     end
    
     % Print on screen the matrices [Cf],[a],[b],[�]
    sp=[Cf',a',b',gt'];
    format short e;
    disp('     [Cf]          [a]          [b]          [�]');        
    disp(sp);
       
    %% Assembling the matrices  A,B,C,D
   
    for i=1:cs
     for j=1:cs 
        
           Ai(i,i)=-1;
           Bi(i,i)=-1; 
           Ci(i,i)=exp(a(i)*Dy(i)); 
           Di(i,i)=-exp(-a(i)*Dy(i));
         if i==j+1
            Ai(i,j)=exp(a(j)*Dy(j));
            Bi(i,j)=exp(-a(j)*Dy(j));
         end;
         if i==j-1 
            Ci(i,j)=-H(i+1)*a(i+1)/a(i);
            Di(i,j)=H(i+1)*a(i+1)/a(i);
         end;
         end
     end
 
    % Corrections added to the matrices
        Ai(1,1)=1;
        Bi(1,1)=1;
        Di(cs,cs)=-Di(cs,cs);
        
    %% Assembling the Final matrix
 
    % Assembling the [�] matrix
        for i=2:cs
        Gt1(i)=gt(i)-gt(i-1); 
        end         
       Gt1(1)=-gt(1);         
       Gt2(cs)=-gt(cs);
 
 
      FM=[Ai,Bi,Gt1';Ci,Di,Gt2'];
         disp('                 ');
         disp('                 ');
         disp('  Final Matrix assembled');
         disp('                 ');
         disp(FM);
                   
    %% Gauss-Jordan Elimination
         GJ=rref(FM);
                 
      rGJ=cs*2;
       cGJ=cs*2+1;
        GJn=GJ(1:rGJ,cGJ);
          format short e;
         
          Ca=GJn(1:cs);
          Cb=GJn(cs+1:rGJ);
            CC=[Ca,Cb];
               disp('                      ') 
               disp('       Ca           Cb');
               disp(CC);
    %% Deferential equation
    Ui=zeros(cs,Sum);
    Qan=0;
    q=0;
    ss=1; 
    for z=1:cs
    
         W=zeros(1,Sum);
             yi=linspace(0,Dy(z),Sum);
                if z==1         
                   y(ss:z*Sum)=yi;
                else
                   y(ss:z*Sum)=yi+y((z-1)*Sum);
                end;
 
     W=Ca(z)*exp(a(z)*yi)+ Cb(z)*exp(-a(z)*yi)+gt(z); 
         Ui(z,1:Sum)=sqrt(W);
              U(ss:z*Sum)=Ui(z,1:Sum);
                  ss=ss+Sum;
                     
    % Calculate the Discharge Analyticaly
        dy=Dy(z)/Sum;
        for l=1:Sum
           q=Ui(z,l)*h(z)*dy;
           Qan=Qan+q;
           q=0;
        end  
        dy=0; 
      end
      
      % $Example
      % Discharge for compound channel with two rectangular cross-sections
      % Manning equation
      A1=h(1)*Dy(1);
      A2=h(2)*Dy(2);
      P1=h(1)+Dy(1);  
      P2=h(2)+Dy(2)+(h(2)-h(1));     
      R1=A1/P1;
      R2=A2/P2;
      
      q1=1/n(1)*A1*R1^(2/3)*s^0.5;
      q2=1/n(2)*A2*R2^(2/3)*s^0.5;
      Qm=q2+q1;
      
      Q=[Qan,Qm];
      
      % Print on screen the results 
      disp('      Discharge (m^3/s)');
      disp(' ');
      disp('       Qan        Qman  ');
      disp(Q);
      % $End of Exmaple
      
      %% Plot 
      B=sum(Dy);
      
      plot(y,U);
      %Setting the graph borders
      axis ([0 B+B/20 0 B+B/20]) ;
      xlabel('y(m)');
      ylabel('U(m/s)');
      title('Compound open-channel velocity profile');
      
      %________________________________## END ##_________________________________
      
