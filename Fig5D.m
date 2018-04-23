function Fig5D()

[T,Y,c] = growth_sim5D(0.1,1);
[T2,Y2,~] = growth_sim5D(0.1,100);

TT = 0:1:3000;
YY = interp1(T,Y(:,11),TT);

figure();
hold on
[AX,H1,H2] = plotyy(TT,YY,0:3000,c(0:3000));
plot(T2,Y2(:,11), 'k--', 'Linewidth',2)

H1.Color = 'k';
H2.Color  = 'r';
H1.LineWidth = 2;
H2.LineWidth = 2;

set(AX(1),'YColor','k','YLim',[0,10],'XLim',[0,3000],'YTick',[0 5 10],'FontSize',20)
set(AX(2),'YColor','r','YLim',[0,0.3],'XLim',[0,3000],'YTick',[0 0.1 0.2 0.3],'FontSize',20)
AX(2).YLabel.String = 'NGF (nM)';
AX(2).YLabel.Interpreter = 'latex';
AX(2).YLabel.FontSize = 24;
AX(1).YLabel.String = '$G^*$ (nM)';
AX(1).YLabel.Interpreter = 'latex';
AX(1).YLabel.FontSize = 24;
AX(1).XLabel.String = 'time (min)';
AX(1).XLabel.Interpreter = 'latex';
AX(1).XLabel.FontSize = 24;


function [T,Y,c] = growth_sim5D(c0,transport_rate)

c = @(t,x) c0 + c0*heaviside(t-500) - c0*heaviside(t-1000);        
cb = c(0,0);

%%% activation / inhibition %%%
At = 1000;
It = 1000;
kpA = 1; kmA = 0.1; kA = kmA/kpA;
kpI = 0.01; kmI = 0.1; kI = kmI/kpI;
kr = transport_rate*0.0025; kdegA = transport_rate*0.005;
ka = transport_rate*0.0025; kdegI = transport_rate*0.005;

%%%% grad sense parameters %%%%
a2 = 0.01; b2 = 1;
a3 = 0.01; %b3 = 0.01;
a4 = a2/It; b4 = 0.00001;
a5 = a3/At*kdegA/kr; b5 = 0.00001;
b3 = b2*kA/kI*a5/a4;
g1 = 0.01; g2 = 0.01;

Ft = 10;
a6 = 1; b6 = 1;
kX = 1; kY = 1;

%%%% growth parameters %%%%
a0 = 0.1; b0 = 0.01;
a1 = 0.000005; b1 = b0/40;
kext = 0.02;
aF = 0.06;

%%% steady-state %%%

A = @(s) At*s./(s+kA);
I = @(s) It*s./(s+kI);
Gstar = @(s1,s2) a1*kdegI/(ka*b1*b0)*(A(s1)./(b0*kdegI/(b1*ka)+(a1*kdegI/(b1*ka))*A(s1)+I(s2)));
Xbar = @(s1,s2) 1/2*(a4*kr/(b2*kdegA)*A(s1).*(a2/a4-I(s2))-b4/g1) + 1/2*sqrt((a4*kr/(b2*kdegA)*A(s1).*(a2/a4-I(s2))-b4/g1).^2 + 4*a2*b4*kr/(b2*kdegA*g1)*A(s1));
Ybar = @(s1,s2) 1/2*(a5/b3*I(s2).*(a3/a5-kr/kdegA*A(s1))-b5/g2) + 1/2*sqrt((a5/b3*I(s2).*(a3/a5-kr/kdegA*A(s1))-b5/g2).^2 + 4*a3*b5/(b3*g2)*I(s2));
ZXbar = @(s1,s2) a4*kr/kdegA*A(s1).*I(s2)./(b4 + g1*Xbar(s1,s2));
ZYbar = @(s1,s2) a5*kr/kdegA*A(s1).*I(s2)./(b5 + g2*Ybar(s1,s2));
J1 = kX/Ft;
J2 = kY/Ft;
v1 = @(s1,s2) a6*Xbar(s1,s2);
v2 = @(s1,s2) b6*Ybar(s1,s2);
B = @ (s1,s2) v2(s1,s2)-v1(s1,s2)+J1*v2(s1,s2)+J2*v1(s1,s2);
Fstar = @(s1,s2) Ft*2*v1(s1,s2).*J2./(B(s1,s2)+sqrt(B(s1,s2).^2-4*(v2(s1,s2)-v1(s1,s2)).*v1(s1,s2)*J2));

FG = @(s1,s2) (a0+aF*Fstar(s1,s2)).*Gstar(s1,s2);

%y = [A,Ac,I,Ic,X,Y,Zx,Zy,Fstar,G,Gstar,R]
f = @(t,y) [kpA*c(t,y(12)).*(At-y(1)) - kmA*y(1);
            kr*y(1) - kdegA*y(2);
            kpI*c(t,0).*(It-y(3)) - kmI*y(3);
            ka*y(3) - kdegI*y(4);
            a2*y(2) - b2*y(5) - g1*y(7).*y(5);
            a3*y(3) - b3*y(6) - g2*y(8).*y(6);
            a4*y(2).*y(3) - b4*y(7) - g1*y(7).*y(5);
            a5*y(2).*y(3) - b5*y(8) - g2*y(8).*y(6);
            a6*y(5).*(Ft - y(9))./(kX + (Ft - y(9))) - b6*y(6).*y(9)./(kY + y(9));
            a0 + aF*y(9) - a1*y(1).*y(10) + b1*y(4).*y(11) - b0*y(10);
            a1*y(1).*y(10) - b1*y(4).*y(11) - b0*y(11);
            0.05+kext*y(11)];

        
IC = [A(cb),kr/kdegA*A(cb),I(cb),ka/kdegI*I(cb),Xbar(cb,cb),Ybar(cb,cb),...
      a4*kr/kdegA*A(cb)*I(cb)./(b4+g1*Xbar(cb,cb)),a5*kr/kdegA*A(cb)*I(cb)./(b5+g2*Ybar(cb,cb)),...
      Fstar(cb,cb),FG(cb,cb)*(b0+b1*ka/kdegI*I(cb))/(a1*A(cb)),FG(cb,cb),0];
     
[T,Y] = ode15s(f,[0,3000],IC);

