x1 = 0;
x2 = 0;
u = 0;
k1 = 9000000;
k2 = 7500;
uin = 20;
uref = 15;
R = 15;
L = 0.0015 ;
C = 0.0001;
iref = (uref*uref+uin*uref)/(uin*R);
pi = 3.14;
n = 1;
t = 0;
Dt = 0.00002;%frequency = 50khz

for i = 1:10000

    ym = L*iref^2/2+C*uref^2/2+C*uin*uref;
    dym = 0;
    ddym = 0;
    
    Dx1 = -x2/L +((uin+x2)/L)*u;
    Dx2 = x1/C -x2/(R*C)-x1*u/C;
    
    x1 = x1+Dx1*Dt;
    x1_store(:,n)=[x1,iref];
    x2 = x2+Dx2*Dt;
    x2_store(:,n)=[x2,uref];
    y = L*x1^2/2+C*x2^2/2+C*uin*x2;
    Dy = L*x1*Dx1+C*x2*Dx2+C*uin*Dx2;
    
    Lf2h = -uin*x2/L-(2*x1*x2+uin*x1)/(C*R)+(2*x2^2+uin*x2)/(R*R*C);
    LgLfh = (uin*uin+x2*uin)/L+(2*x1*x2+uin*x1)/(R*C);
   v = ddym - k2*(Dy - dym)-k1*(y-ym);

   u = (v-Lf2h)/LgLfh;
    
    y_store(:,n)=[y;ym];

   t=t+Dt;
    n=n+1;  
end

figure(1)
plot((1:n-1)*Dt,y_store)
xlabel('time/s')
ylabel('y')
legend('y-act','y-ref')
figure(2)
plot((1:n-1)*Dt,x1_store)
xlabel('time/s')
ylabel('x1')
legend('x1-act','x1-ref')
figure(3)
plot((1:n-1)*Dt,x2_store)
xlabel('time/s')
ylabel('x2')
legend('x2-act','x2-ref')













