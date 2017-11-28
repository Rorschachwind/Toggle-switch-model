%{
11/21/2017 latest version of code
Toggle switch genetic model
%}

clear;clc;
global d kmax1 K2 Ot3 K5 kmax2 K3 k0 d1 
%degradation rate constant
d = 1.0;
d1 = 0.3;
%max rate
kmax1 = 200.0;
kmax2 = 10.0;
% operator concentration
Ot3 = 1.0;
%basal level synthetic rate
k0 = 0.01;
% hill, synthesis constant
K = 0.5:0.5:3;
% K2 = 1.0;
K5 = 1.0;
K3=2.0;

%% test y^2=F(x)
% x nullcline, y^2 = F(x)
% first draw F(x) to see where it is positive
xmin = k0/d;	% for x<xmin, F<0, at x=xmin, F has an asymptote
% 		% so, no need to look at negative x
x = [-2:0.001:0.9*xmin 1.01*xmin:0.001:2];
K2 = 1.0;
A = kmax1*K2*Ot3;	% combine three constants
F = (1/K5)*( (A*x.^2)./((1+K2*(x.^2)).*(d*x-k0)) - 1);
figure(1); 
plot(x,F,'.-'); grid on;
xlabel('X-NRI','Fontsize',15);ylabel('Y-LacI','Fontsize',15);
title('F(x)=y^2 changes with X','Fontsize',15);

%% Construct phase plane and time courses for different K2
tspan = [0 3000.0];
% y0 = [0.52; 11.7; ]; 5this initial condition is closer to fixed point
y0=[2.0; 1.0;];
x0=y0(1);
% K2 vector, which can be changed to different region
KKK=3:1:8;
n=length(KKK);
for i=1:n
    K2=KKK(i);
    [t,y] = ode45(@repress_fun, tspan, y0);
    
% time courses at different K2
     figure(2);
     subplot(3,2,i);
     plot(t,y(:,1),'b',t,y(:,2),'r');
     ylim([0 20]);
     title(['K2= ', num2str(K2)]);

%add fixed point
    A= kmax1*K2*Ot3;
    B=kmax2*K3;
    
%new method to compute intersection
    x=[0:0.001:0.9*xmin 1.01*xmin:0.001:1];
    y1=d*x-k0;
    y2=d1^2*A*x.^2.*(1+K3*x.^2).^2./((K5*B^2*x.^4+d1^2*(1+K3*x.^2).^2).*(1+K2*x.^2));
    %use a new function InterX
    rts=InterX([x;y1],[x;y2]);
    rts=rts(1,:);
    Y=(kmax2*K3.*rts.^2)./(d1*(1+K3*rts.^2));
    m=length(rts);   
    for j=1:m
             figure(i+2);
             plot(rts(j),Y(j),'o');
             hold on;
    end
%add nullcline
%{   
This vector is not right, as it will take k0/d into x
x1=-5;x2=80;dx=(x2-x1)/10;
y1=-10;y2=80;dy=(y2-y1)/10;

The vector used below has excluded x=k0/d
%}
     figure(i+2);
     x1=1.01*xmin;x2=15;dx=(x2-x1)/10; %here only plot when x >xmin, then xnullcline can be solved
     y1=0;y2=35;dy=(y2-y1)/10;
     vx=x1:dx/10000:x2;
     vy=y1:dy:y2;
     vxx=x1:dx:x2;
     ynull= (kmax2*K3.*vx.^2)./(d1*(1+K3*vx.^2));
     xnull= ((K2*kmax1*Ot3*vx.^2)./(K5*(d*vx-k0).*(1+K2*vx.^2))-1/K5).^(1/2);
     plot(vx,xnull,'y--',vx,ynull,'g--'); grid on; hold on;

%add vector field
     figure(i+2);
     [U,V]=meshgrid(vxx,vy);
     V_x=k0+(K2*kmax1*Ot3*U.^2)./((1+K2*U.^2).*(1+K5*V.^2))-d*U;
     V_y=(kmax2*K3*U.^2)./(1+K3*U.^2)-d1*V;
     quiver(U,V,V_x,V_y);
% add trajectory
     plot(y(:,1),y(:,2));
     hold off;
     axis([x1 x2 y1 y2]);
     xlabel('X-NRI','Fontsize',15);ylabel('Y-LacI','Fontsize',15);
     title(['Phase plane of system','K2= ',num2str(K2)],'Fontsize',18);
     legend('fixed point','xnull','ynull','vector field','trajectory');


end

figure(2);
xlabel('Time/min');ylabel('Concentration');
legend('NRI','LacI');
suptitle('Time course for NRI, LacI at different K2');


%% delta_tau graph
KK = 3:1:10;
figure(9);
for i=1:length(KK)
    K2=KK(i);
    A= kmax1*K2*Ot3;
    B=kmax2*K3;
    
    x=[0:0.001:0.9*xmin 1.01*xmin:0.001:1];
    y1=d*x-k0;
    y2=d1^2*A*x.^2.*(1+K3*x.^2).^2./((K5*B^2*x.^4+d1^2*(1+K3*x.^2).^2).*(1+K2*x.^2));
%use a new function InterX
    rts=InterX([x;y1],[x;y2]);
    rts=rts(1,:);
    Y=(kmax2*K3.*rts.^2)./(d1*(1+K3*rts.^2));
    m=length(rts)
    kk=0:1e-3:5;
    t1=sqrt(4*kk);
    
% determine type of fixed pts.
    for j=1:m
        if ((isreal(rts(j))&&(rts(j)>=0)))
            h=rts(j);
            hh=Y(j);
            fx= -d+(2*A*h)/((1+K5*hh^2)*(1+K2*h^2)^2);
            fy= -2*A*K5*h^2*hh/((1+K2*h^2)*(1+K5*hh^2)^2);
            gx= 2*B*h/((1+K3*h^2)^2);
            gy= -d1;
            AA=[fx fy;gx gy];
            eig(AA)
            del=det(AA);
            tau=trace(AA);
            
            if (del<0)
                  disp(['K2= ', num2str(K2),'fixed point= ',num2str(h),'type= ','unstable saddle']);
            elseif ((tau>0)&&(del<tau^2/4))
                disp(['K2= ', num2str(K2),'fixed point= ',num2str(h),'type= ','unstable node']);
            elseif ((tau>0)&&(del>tau^2/4))
                disp(['K2= ', num2str(K2),'fixed point= ',num2str(h),'type= ','unstable spiral']);
            elseif ((tau<0)&&(del<tau^2/4))
                disp(['K2= ', num2str(K2),'fixed point= ',num2str(h),'type= ','stable node']);
            else
                disp(['K2= ', num2str(K2),'fixed point= ',num2str(h),'type= ','stable spiral']);
            end
           
        end
    end

% plot three del_tau plane
    if m ==3
        for j=1:m
            h=rts(j);
            hh=Y(j);
            fx= -d+(2*A*h)/((1+K5*hh^2)*(1+K2*h^2)^2);
            fy= -2*A*K5*h^2*hh/((1+K2*h^2)*(1+K5*hh^2)^2);
            gx= 2*B*h/((1+K3*h^2)^2);
            gy= -d1;
            AA=[fx fy;gx gy];
            delt=det(AA);
            tau=trace(AA);

            subplot(2,2,j);
            plot(kk,t1,'r--',kk,-t1,'r--',delt,tau,'k.');grid on;hold on;
            axis([-1 3 -4 4]); 
            xlabel('delta'); ylabel('tau');
            title(['fixed point ' num2str(j)]);
        end
    elseif m == 1
            h=rts(1);
            hh=Y(1);
            fx= -d+(2*A*h)/((1+K5*hh^2)*(1+K2*h^2)^2);
            fy= -2*A*K5*h^2*hh/((1+K2*h^2)*(1+K5*hh^2)^2);
            gx= 2*B*h/((1+K3*h^2)^2);
            gy= -d1;
            AA=[fx fy;gx gy];
            delt=det(AA);
            tau=trace(AA);
        if rts(1)<=0.05
            subplot(2,2,1);
            plot(kk,t1,'r--',kk,-t1,'r--',delt,tau,'k.');grid on;hold on;
        end
        if rts(1)>0.1
            subplot(2,2,3);
            plot(kk,t1,'r--',kk,-t1,'r--',delt,tau,'k.');grid on;hold on;
        end
    end

 end

subplot(2,2,1);
axis([-1 3 -4 4]); 
xlabel('delta','Fontsize',14); ylabel('tau','Fontsize',14);
title('fixed point 1','Fontsize',14);
subplot(2,2,2);
axis([-1 3 -4 4]); 
xlabel('delta','Fontsize',14); ylabel('tau','Fontsize',14);
title('fixed point 2 ','Fontsize',14);
subplot(2,2,3);
axis([-1 3 -4 4]); 
xlabel('delta','Fontsize',14); ylabel('tau','Fontsize',14);
title('fixed point 3','Fontsize',14 );
%% bifurcation diagram
% set region of K2 to acquire enough point
KK2 = 0.01:0.0005:1;
m=length(KK2);
figure(10);
for i=1:m
    K2=KK2(i);
    A= kmax1*K2*Ot3;
    B=kmax2*K3;
    x=[0:0.001:0.9*xmin 1.01*xmin:0.001:1];
    y1=d*x-k0;
    y2=d1^2*A*x.^2.*(1+K3*x.^2).^2./((K5*B^2*x.^4+d1^2*(1+K3*x.^2).^2).*(1+K2*x.^2));
    %use a new function InterX
    rts=InterX([x;y1],[x;y2]);
    rts=rts(1,:);
    Y=(kmax2*K3.*rts.^2)./(d1*(1+K3*rts.^2));
    m=length(rts);   
    for j=1:m
             plot(K2,rts(j),'k.');
             hold on;
    end
      
end
xlabel('K2','Fontsize',15);ylabel('X coordinate of fixed point','Fontsize',15);
title('Bifurcation Diagram','Fontsize',15);
