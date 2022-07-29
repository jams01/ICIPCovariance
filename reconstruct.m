load('ready/H')
load('ready/tf')
load('ready/P')
load('ready/params')

%load measurements
load('medidas/Cuento_aperture8rand=false.mat')

% hand crafted correction
X=imrotate(X,-0.6,'bilinear');
% select the same region chosen in the sensing matrix 
sH = X(513:513+params.M-1,64:64+params.N-1,:);
Y = zeros([round(M/disminuir),(N/disminuir),shots]); 
% decimate
for si=1:shots
    Y(:,:,si) = imresize(sH(:,:,nshots(si)),[round(M/disminuir),(N/disminuir)],'bilinear');
end

Y1 = reshape(Y,[round(M/disminuir)*round(N/disminuir),shots]);
% extract the partition of the image using tf (the spatial distribution of the filters)
tf = reshape(tf,[round(M/disminuir)*round(N/disminuir),1]);
Ys = cell(1,patterns);
for i=1:patterns
    Ys{i}=Y1(tf==i,:)';
end

for i=1:patterns
    P{i}=P{i}';
end
% we compute the mean and substract it from measurements, and occording
% repeat this process since it gives better results.
meanst=estimate_mean(P,Ys);
Y1=cell(size(Ys));
for i=1:patterns
    Y1{i}=Ys{i}-P{i}'*(kron(meanst,ones(1,size(Ys{i},2))));
end
for j=1:15
    mst=estimate_mean(P,Y1);
    meanst=meanst+mst;
    for i=1:patterns
        Y1{i}=Y1{i}-P{i}'*(kron(mst,ones(1,size(Ys{i},2))));
    end
end
it = 100;
%estimate the covariance matrix
Sigma1 =  estimate_cov6(P,Y1,it,0.0005,200,4,0,0,0);
% extract the eigenvectors
[Wrm,e1]=svd(Sigma1);
Wrm=real(Wrm);
% select only 5
W2=Wrm(:,1:5);
Y_est = cell(1,patterns);
% reconstruct
for i=1:patterns
    Y_est{i}=pinv(P{i}'*W2)*Y1{i};%
end
Y_est=W2*cell2mat(Y_est);%
Ys1 = zeros(floor(101/disminuirbands), M/disminuir*N/disminuir);
n = 1;
% reorder
for i=1:patterns
    Ys1(:,tf==i)=(Y_est(:,n:((n-1)+sum(tf==i))));
    n = n + sum(tf==i);
end
Ys1=Ys1+meanst;
imhyp = reshape(Ys1',[M/disminuir,N/disminuir,floor(101/disminuirbands)]);

% the reconstructed hyperspectral image
imhyp=max(0,imhyp);




% a signature of the white to white balance
load('white_reference.mat')
b1=reshape(meanst1,[1,1,floor(101/disminuirbands)]);
b1=b1./max(b1);
imhyp1=imhyp;
imhyp1 = imhyp1 ./ b1;

%create a integration curve to emulate an rgb camera
integration = normpdf([-2:.2:2.1],0,1);
integration=integration./sum(integration);
%integration = sqrt(integration);
integration(:) = integration / sum(integration);
imt =zeros(M/disminuir,N/disminuir,3);
    
    int1 = integration(4:17)./sum(integration(4:17));
    int1(1:3)=0;
    imt(:,:,3) =  imgaussfilt(sum(imhyp1(:,:,1:14).*reshape(int1,[1,1,14]),3),1);
    int1 = integration(6:15)./sum(integration(6:15));
    imt(:,:,2) =  imgaussfilt(sum(imhyp1(:,:,15:24).*reshape(int1,[1,1,10]),3),1);
    int1 = integration(6:12)./sum(integration(6:12));
    imt(:,:,1) =  imgaussfilt(sum(imhyp1(:,:,25:end).*reshape(int1,[1,1,7]),3),1);
    
figure,imagesc(imt/12)
