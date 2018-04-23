function Fig5AB()

%% data %%
load('./grad_data.mat')
cf = [0.92, 1.1, 1.5, 2.1];
cn = [0.001, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100];
s = [0.12 , 0.18, 0.24, 0.3];


outgrowth{length(s), length(cn)} = [];
bias{length(s), length(cn)} = [];
n_val_og(length(s), length(cn)) = 0;
n_val_bs(length(s), length(cn)) = 0;

for j= 1:length(s)
    for k = 1:length(cn)
        jk = find((grad_data.NGF == cn(k)) & (grad_data.gradient == s(j)));
        outgrowth{j,k} = grad_data.averageOutgrowth(jk);
        bias{j,k} = grad_data.directionalBias(jk);
        bias{j,k}(outgrowth{j,k} < 100)=[];
        n_val_og(j,k) = length(outgrowth{j,k});
        n_val_bs(j,k) = length(bias{j,k});
    end
end

og_mean = cell2mat(cellfun(@mean, outgrowth, 'uni', 0));
og_std = cell2mat(cellfun(@std, outgrowth, 'uni', 0));
bs_mean = cell2mat(cellfun(@mean, bias, 'uni', 0));
bs_std = cell2mat(cellfun(@std, bias, 'uni', 0));


%% model %%
s = [0.12 , 0.18, 0.24, 0.3]*1e-3;
cc = logspace(-4,3,100);
lin = {'b-','r-','g-','k-'};
mar = {'bo','ro','go','ko'};

for ind = 1:4

    for k = 1:length(cc)
        [OG(k), BS(k)] = simplant(cc(k),s(ind));
    end

    figure(1)
    hold on
    
    g{1} = plot(cc, OG, lin{ind}, 'LineWidth', 2);
    errorbar(cn.*cf(ind), og_mean(ind,:), og_std(ind,:)./sqrt(n_val_og(ind,:)),...
        mar{ind}, 'MarkerSize', 8, 'LineWidth', 1.5)
    
    xlabel('NGF (nM)', 'Interpreter', 'Latex', 'FontSize', 22)
    ylabel('radial outgrowth ($\mu{m}$)', 'Interpreter', 'Latex', 'FontSize', 22)
    set(gca,'XScale', 'log', 'FontSize', 20,'xtick',[0.001,0.01,0.1,1,10,100],'xticklabel',[0.001,0.01,0.1,1,10,100])
    axis([0.0002,500,0,1000])
    
    figure(2)
    hold on
    
    h(ind)=plot(cc , BS, lin{ind}, 'LineWidth', 2);
    errorbar(cn.*cf(ind), bs_mean(ind,:), bs_std(ind,:)./sqrt(n_val_bs(ind,:)),...
        mar{ind},'MarkerSize',8,'LineWidth',1.5)
    
    xlabel('NGF (nM)', 'Interpreter', 'Latex', 'FontSize', 22)
    ylabel('directional bias', 'Interpreter', 'Latex', 'FontSize', 22)
    set(gca,'XScale', 'log', 'FontSize', 20,'xtick',[0.001,0.01,0.1,1,10,100],'xticklabel',[0.001,0.01,0.1,1,10,100])
    axis([0.0002,500,-0.1,0.2])

end
legend(h,{'0.12 %','0.18 %','0.24 %','0.3 %'})

%% sub-functions %%

function [outgrowth,bias] = simplant(c0,s)

tt = linspace(0,2*pi,361);
tt(end)=[];
R = zeros(1,length(tt));

for k=1:length(tt)
    R(k) = growth_sim5A(tt(k),c0,s);
end

[A, B] = coeffs(R);
outgrowth = A(1);
bias = B(2)/A(1);

function R = growth_sim5A(th,c0,s)

rE = 300;
c = @(t,x) c0*exp(s*(rE + x)*sin(th));
cb = c(0,0);

%%% activation / inhibition %%%
At = 1000;
It = 1000;
kpA = 1; kmA = 0.1; kA = kmA/kpA;
kpI = 0.01; kmI = 0.1; kI = kmI/kpI;
kr = 0.0025; kdegA = 0.005;
ka = 0.0025; kdegI = 0.005;

%%%% grad sense parameters %%%%
a2 = 0.01; b2 = 1;
a3 = 0.01; %b3 = 0.01;
a4 = a2/It; b4 = 0.00001;
a5 = a3/At*kdegA/kr; b5 = 0.00001;
b3 = b2*kA/kI*a5/a4;
g1 = 0.01; g2 = 0.01;

%%%% Golbeter-Koshland %%%%
Ft = 10;
a6 = 1; b6 = 1;
kX = 1; kY = 1;

%%%% growth parameters %%%%
a0 = 0.1; b0 = 0.01;
a1 = 0.000005; b1 = b0/40;
kext = 0.02; r0 = 0.05;
aF = 0.06;

%%% steady-state %%%
A = @(s) At*s./(s+kA);
I = @(s) It*s./(s+kI);
Gstar = @(s1,s2) a1*kdegI/(ka*b1*b0)*(A(s1)./(b0*kdegI/(b1*ka)+(a1*kdegI/(b1*ka))*A(s1)+I(s2)));
Xbar = @(s1,s2) 1/2*(a4*kr/(b2*kdegA)*A(s1).*(a2/a4-I(s2))-b4/g1) +...
    1/2*sqrt((a4*kr/(b2*kdegA)*A(s1).*(a2/a4-I(s2))-b4/g1).^2 + 4*a2*b4*kr/(b2*kdegA*g1)*A(s1));
Ybar = @(s1,s2) 1/2*(a5/b3*I(s2).*(a3/a5-kr/kdegA*A(s1))-b5/g2) +...
    1/2*sqrt((a5/b3*I(s2).*(a3/a5-kr/kdegA*A(s1))-b5/g2).^2 + 4*a3*b5/(b3*g2)*I(s2));
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
            r0+kext*y(11)];

        
IC = [A(cb),kr/kdegA*A(cb),I(cb),ka/kdegI*I(cb),Xbar(cb,cb),Ybar(cb,cb),...
      a4*kr/kdegA*A(cb)*I(cb)./(b4+g1*Xbar(cb,cb)),a5*kr/kdegA*A(cb)*I(cb)./(b5+g2*Ybar(cb,cb)),...
      Fstar(cb,cb),FG(cb,cb)*(b0+b1*ka/kdegI*I(cb))/(a1*A(cb)),FG(cb,cb),0];
%     
[~,Y] = ode15s(f,[0,2880],IC);
R = Y(end,12);

function [A,B] = coeffs(X)

tt = linspace(0,2*pi,361);
tt(end) = [];

Y = zeros(1,length(tt));
N = length(tt);
for k = 1:length(Y)
    b = exp(-2*pi*1i/N*(k-1)*(0:N-1));
    Y(k) = b*X';                
end    
Y = Y/N;
Y(2:180)=2*Y(2:180);
Y(181:360)=[];

A = real(Y);
B = -imag(Y);