function Fig3D()

noise = 0.1;
trials = 100;
rng(0)

At = 1000;
It = 1000;

kA = 0.1;
kI = 10;

A = @(s) At*s./(s+kA);
I = @(s) It*s./(s+kI);

ss = logspace(-3,2,100);

for kk = 1:trials
    for ii = 1:100
        for jj = 1:100
            
    %%%% gradient sensor %%%%%%
    a2 = 0.01; b2 = 1;
    a3 = 0.01;
    a4 = a2/It; b4 = 0.00001;
    a5 = a3/At; b5 = 0.00001;
    b3 = b2*kA/kI*a5/a4;
    g1 = 0.01; g2 = 0.01;

    %perturbations to parameters
    a2 = a2*(1-noise+rand*(2*noise));
    b2 = b2*(1-noise+rand*(2*noise));
    a3 = a3*(1-noise+rand*(2*noise));
    b3 = b3*(1-noise+rand*(2*noise));
    
    Xbar = @(s1,s2) 1/2*(a4/b2*A(s1).*(a2/a4-I(s2))-b4/g1) +...
        1/2*sqrt((a4/b2*A(s1).*(a2/a4-I(s2))-b4/g1).^2 + 4*a2*b4/(b2*g1)*A(s1));
    Ybar = @(s1,s2) 1/2*(a5/b3*I(s2).*(a3/a5-A(s1))-b5/g2) +...
        1/2*sqrt((a5/b3*I(s2).*(a3/a5-A(s1))-b5/g2).^2 + 4*a3*b5/(b3*g2)*I(s2));
    ZXbar = @(s1,s2) a4*A(s1).*I(s2)./(b4 + g1*Xbar(s1,s2));
    ZYbar = @(s1,s2) a5*A(s1).*I(s2)./(b5 + g2*Ybar(s1,s2));

    %goldbetter-koshland parameters
    Ft = 10;
    a6 = 1; b6 = 1;
    Km1 = 1; Km2 = 1;
    J1 = Km1/Ft; J2 = Km2/Ft;
    
    v1 = @(s1,s2) a6*Xbar(s1,s2);
    v2 = @(s1,s2) b6*Ybar(s1,s2);
    B = @(s1,s2) v2(s1,s2)-v1(s1,s2)+J1*v2(s1,s2)+J2*v1(s1,s2);
    Fstar = @(s1,s2) 2*v1(s1,s2).*J2./(B(s1,s2)+sqrt(B(s1,s2).^2-4*(v2(s1,s2)-v1(s1,s2)).*v1(s1,s2)*J2));

    Z(kk,ii,jj) = (Fstar(ss(ii),ss(jj)));
        end
    end
end


%% plotting %%
Z = squeeze(mean(Z,1));
figure()
imagesc(Z);
set(gca, 'ydir', 'normal');
cbh = colorbar;
set(cbh, 'YTick', [0:0.25:1])
ylabel(cbh, '$\bar{F^*}(A(c_1),I(c_2))$', 'Interpreter', 'Latex', 'FontSize', 20)
caxis([0 1])

ticks = [1,20,40,60,80,100];
labels = {'0.001','0.01','0.1','1','10','100'};
set(gca, 'XTick', ticks, 'XTickLabel', labels);
set(gca, 'YTick', ticks, 'YTickLabel', labels);
xlabel('$c_2$ (nM)','Interpreter','Latex','FontSize',20)
ylabel('$c_1$ (nM)','Interpreter','Latex','FontSize',20)
set(gca,'FontSize',16)
axis square
colormap gray

