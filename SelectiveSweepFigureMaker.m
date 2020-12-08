
%PG, Ps, Ng, Ns

c=0.01;
epsilon=10^-2;
fa=1.01;
fb=1;

conjRate = @(y) c*y(1)+epsilon*c*y(2);
fbar = @(y) [fa,fa,fb,fb]*y;

Deriv= @(t,y) diag([fa,fa,fb,fb]-fbar(y))*y + [0,0,conjRate(y),0;0,0,0,conjRate(y);0,0,-conjRate(y),0;0,0,0,-conjRate(y) ]*y;

tspan = [0 1000];
y0 = [10^-3;0;0.9-10^-3;0.1];
y0=y0./sum(y0);

cw= c*(y0(3) +epsilon*y0(4))./(y0(3) +y0(4));
cm= c*(y0(1) +epsilon*y0(2))./(y0(1) +y0(2));
del_f=fa-fb;
M0=y0(1)+y0(2);
Minf=(cm/cw).^(-del_f/(del_f+cw))* (M0^(cw/(del_f+cw))) ;

Ginf= Minf*(y0(1))./(y0(1) +y0(2)) + (1-Minf)*(y0(3))./(y0(3) +y0(4));
Sinf=1-Ginf;

subplot(1,2,1);
[t,y] = ode45(Deriv, tspan, y0);
plot(t,y(:,1), 'Color',[0.3,0.3,1],'lineWidth',2);
hold on
plot(t,y(:,2),':', 'Color',[0.2,1,0.2],'lineWidth',2);
plot(t,y(:,3), 'Color',[0,0,0.5],'lineWidth',2);
plot(t,y(:,4),':','Color',[0,0.5,0],'lineWidth',2);

scatter(t(end),Ginf,'dk');
scatter(t(end),Sinf,'dk');


k= -(cm/cw)* (M0)^(-cw/del_f);

alpha=(del_f+cw)/del_f;
%%Newton Solve!
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;

Ginf= Minf*(y0(1))./(y0(1) +y0(2)) + (1-Minf)*(y0(3))./(y0(3) +y0(4));
Sinf=1-Ginf;
scatter(t(end),Ginf,'ko');
scatter(t(end),Sinf,'ko');
%PG, Ps, Ng, Ns
legend({'G_P','S_P','G_N','S_N'})

ylim([0,1])



y0 = [0;10^-3;0.9-10^-3;0.1];
y0=y0./sum(y0);

subplot(1,2,2);
[t,y] = ode45(Deriv, tspan, y0);
plot(t,y(:,1), 'Color',[0.3,0.3,1],'lineWidth',2);
hold on
plot(t,y(:,2),':', 'Color',[0.2,1,0.2],'lineWidth',2);
plot(t,y(:,3), 'Color',[0,0,0.5],'lineWidth',2);
plot(t,y(:,4),':','Color',[0,0.5,0],'lineWidth',2);

cw= c*(y0(3) +epsilon*y0(4))./(y0(3) +y0(4));
cm= c*(y0(1) +epsilon*y0(2))./(y0(1) +y0(2));
del_f=fa-fb;
M0=y0(1)+y0(2);
Minf=(cm/cw).^(-del_f/(del_f+cw))* (M0^(cw/(del_f+cw))) ;

Ginf= Minf*(y0(1))./(y0(1) +y0(2)) + (1-Minf)*(y0(3))./(y0(3) +y0(4));
Sinf=1-Ginf;

scatter(t(end),Ginf,'dk');
scatter(t(end),Sinf,'dk');

k= -(cm/cw)* (M0)^(-cw/del_f);

alpha=(del_f+cw)/del_f;
%%Newton Solve!
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;
val= 1+k*Minf^alpha + Minf*(cm-cw)/cw;
deriv= alpha*k*Minf^(alpha-1) +(cm-cw)/cw;
Minf=Minf - val/deriv;

Ginf= Minf*(y0(1))./(y0(1) +y0(2)) + (1-Minf)*(y0(3))./(y0(3) +y0(4));
Sinf=1-Ginf;
scatter(t(end),Ginf,'ko');
scatter(t(end),Sinf,'ko');


%PG, Ps, Ng, Ns
ylim([0,1])
