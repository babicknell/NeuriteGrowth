function S =  explant_param()

%Script for parameterising the boundary curves of binarised images of DRG 
%explants that have been segmented into neurite outgrowth and cell-body 
%regions. Area and length-based measuremts of neurite outgrowth are also 
%computed, and the results saved in the data structure 'S'.

filenames = dir('./neurite_masks/*n*.tif');
filenames = {filenames.name};
numFiles = length(filenames);

scale = 2.58; %microns per pixel image scale factor

con = [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100];
conString = {'000p001','000p003','000p010','000p030','000p100','000p300',...
             '001p000','003p000','010p000','030p000','100p000'};
grad = [0, 0.1, 0.2, 0.3, 0.4];
gradString = {'0p0', '0p1', '0p2', '0p3', '0p4'};

%Pre-allocation of data structure fields
Name{numFiles} = [];
Grad{numFiles} = [];
Con{numFiles} = [];
OG{numFiles} = [];
GR{numFiles} = [];
Neurites{numFiles} = [];
Explant{numFiles} = [];
Neu_rho{numFiles} = [];
Exp_rho{numFiles} = [];
Cen{numFiles} = [];
Neurites_F{numFiles} = [];
Explant_F{numFiles} = [];
avOutgrowth{numFiles} = [];
direcBias{numFiles} = [];

for k = 1:numFiles
     
    gradient = filenames{k}(10:12);
    concentration = filenames{k}(14:20);
    Grad{k} = grad(strcmp(gradient,gradString));
    Con{k} = con(strcmp(concentration,conString));
    
    %Load binary masks of neurite (n) and explant (x) regions.
    In = imread(['./neurite_masks/',filenames{k}]);
    Ix = imread(['./explant_masks/',filenames{k}(1:end-5),'x.tif']);
    
    %Reformat image matrices and fit boundary curves
    BWx = Ix > 0;
    BWx = flip(BWx)';
    BWx([1,end],:) = 0;
    BWx(:,[1,end]) = 0;
    
    BW = In > 0;
    BW = flip(BW)';
    BW([1,end],:) = 0;
    BW(:,[1,end]) = 0;
    BW(BWx==1) = 1;
    
    CCo = bwconncomp(BW);
    numPixelso = cellfun(@numel,CCo.PixelIdxList);
    [~,idxo] = max(numPixelso);
    BWn = false(size(BW));
    BWn(CCo.PixelIdxList{idxo}) = true;

    Bx = bwboundaries(BWx);
    Bn = bwboundaries(BWn);
    lx = cellfun(@length,Bx);
    [~,ix] = sort(lx,'descend');
    ln = cellfun(@length,Bn);
    [~,in] = sort(ln,'descend');
    b_x = Bx{ix(1)};
    b_n = Bn{in(1)};
    
    s = regionprops(BWx,'centroid');
    cen = flip(s.Centroid);
    
    %Compute Outgrowth and Guidance Ratio measures
    Neu = BWn - BWx;
    H = Neu; L = Neu;
    H(:,1:round(cen(2))) = 0;
    L(:,round(cen(2)):end) = 0;
    
    OG{k} = sum(Neu(:))/sum(BWx(:));
    GR{k} = (sum(H(:)) - sum(L(:)))/(sum(H(:)) + sum(L(:)));
    
    %Parameterise boundaries for Fourier analysis
    X_n = (b_n(:,1) - cen(1))*scale;
    Y_n = (b_n(:,2) - cen(2))*scale;
    X_x = (b_x(:,1) - cen(1))*scale;
    Y_x = (b_x(:,2) - cen(2))*scale;

    Neurites{k} = [X_n,Y_n];
    Explant{k} = [X_x,Y_x];
    
    %Smooth neurite outgrowth boundary 
    N = min(length(X_n)-1,300);
    X_n = [X_n(end-N:end);X_n;X_n(1:N)];
    Y_n = [Y_n(end-N:end);Y_n;Y_n(1:N)];
    X_n = smooth(X_n,floor(N/2),'moving');
    Y_n = smooth(Y_n,floor(N/2),'moving');
    X_n = X_n(N+1:end-N);
    Y_n = Y_n(N+1:end-N);

    %Convert to polar paramerterisation about centroid
    tt = linspace(0,2*pi,361);
    tt(end) = [];
    Rho_N = zeros(1,length(tt));
    Rho_X = zeros(1,length(tt));

    [theta_x,rho_x] = cart2pol(X_x,Y_x);
    [theta_n,rho_n] = cart2pol(X_n,Y_n);
    theta_x(theta_x<0) = theta_x(theta_x<0)+2*pi; 
    theta_n(theta_n<0) = theta_n(theta_n<0)+2*pi; 
    [phi_x,i_x]=sort(theta_x,'ascend');
    [phi_n,i_n]=sort(theta_n,'ascend');
    rho_x = rho_x(i_x);
    rho_n = rho_n(i_n);

    %Choose closest point to origin when multiply defined
    [phi_xU] = unique(phi_x);
    [phi_nU] = unique(phi_n);
    rho_xU = zeros(1,length(phi_xU));
    rho_nU = zeros(1,length(phi_nU));  
    for m = 1:length(phi_xU)
        rmin_x = min(rho_x(phi_x==phi_xU(m)));    
        rho_xU(m) = rmin_x;
    end
    for m = 1:length(phi_nU)
        rmin_n = min(rho_n(phi_n==phi_nU(m)));    
        rho_nU(m) = rmin_n;
    end
    
    for j = 1:length(tt)
        ind_x = find(phi_xU >= tt(j),1);
        ind_n = find(phi_nU >= tt(j),1);

        if isempty(ind_x)
            ind_x = length(phi_xU);
        end
        if isempty(ind_n)
            ind_n = length(phi_nU);
        end
        
        Rho_X(j) = rho_xU(ind_x); 
        Rho_N(j) = rho_nU(ind_n); 
    end

    %Fill cells for data structure
    Exp_rho{k} = Rho_X;
    Neu_rho{k} = Rho_N;
    Cen{k} = cen;
    Name{k} = filenames{k}(1:end-6);
     
    %Fourier transforms
    Neurites_F{k} = exdft(Neu_rho{k}-Exp_rho{k});
	Explant_F{k} = exdft(Exp_rho{k});
    
    Neurites_F{k} = [real(Neurites_F{k})', -imag(Neurites_F{k})'];
    Explant_F{k} = [real(Explant_F{k})', -imag(Explant_F{k})'];
    
    avOutgrowth{k} = Neurites_F{k}(1,1);
    direcBias{k} = Neurites_F{k}(2,2)/Neurites_F{k}(1,1);
    
end

S = struct('name',Name,'gradient',Grad,'concentration',Con,'outgrowthBoundary',...
    Neurites,'somataBoundary',Explant,'centroid',Cen,'outgrowthCoeffs',...
    Neurites_F,'somataCoeffs',Explant_F,'averageOutgrowth',avOutgrowth,...
    'directionalBias',direcBias,'OG',OG,'GR',GR);

save('S','S')

function [Y]=exdft(X)
    
tt = linspace(0,2*pi,361);
tt(end) = [];

Y = zeros(1,length(tt));
N = length(tt);
for k = 1:length(Y)
    b = exp(-2*pi*1i/N*(k-1)*(0:N-1));
    Y(k) = b*X';                
end    
Y = Y/N;
Y(2:180) = 2*Y(2:180);
Y(182:360) = [];