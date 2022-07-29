

M = 512; %sensor is bigger, but the scene is well focused in the center of the sensor
N = 512; % so, we extract a patch
disminuir = 2; % decimation factor
patterns = 64; % how many patterns search in the sensing matrix
disminuirbands=3.25; % we captured 101 bands, but the prism dispersion was less, this factor reduces to 31 bands
nshots = [4,13,16]; % we capture 19 shots, this ones were the best calibrated
shots=length(nshots);
H = single(zeros(M/disminuir,M/disminuir,floor(101/disminuirbands),shots));
offset = -4;
for k=1:shots
    load(['codigos/False_8/',num2str(nshots(k)),'_aperture8','rand=false']) % load the sensing matrix (captured in lab)
    [n1,n2,n3] = size(X);
    X=imrotate(X,-0.3,'bilinear'); % correction
    sH = X(513:513+M-1,64+offset:64+N+offset-1,:);% Extract the patch

    G=zeros(round(M/disminuir),(N/disminuir),101);
    for si=1:101
        G(:,:,si) = imresize(sH(:,:,si),[round(M/disminuir),(N/disminuir)],'bilinear');%decimate the data
    end
    if disminuirbands>1
        G1 = zeros(M/disminuir,N/disminuir,floor(101/disminuirbands));
        index=1;
        bandsall =450:2:650;
        bandsdisc = zeros(1,floor(101/disminuirbands));
        load('disper31')% dispersion curve captured at lab
        for i=1:floor(101/disminuirbands) % integrate the spectral bands following the dispersion
            G1(:,:,index) = sum(G(:,:,disper==i),3)./sum(disper==i);
            [v,p]=find(disper==i);
            bandsdisc(i)=mean(bandsall(p));
            index=index+1;
        end
    else
        G1=G;
    end
    clear G;
    % some hand crafted corrections
    G1 = real(G1.^1.2);
    G1(70:100,35:50,:) = repmat(median(G1(70:100,35:50,:),1),[31,1]);
    G1(234:236,70:75,:) = repmat(median(G1(234:236,70:75,:),1),[3,1]);
    G1 = G1./max(G1(:));
    H(:,:,:,k)=single(G1);
end
clear sH X

t1 = H(20,1:patterns,:,1); % take a section of the code to look for repeated codes
t1 = normc(reshape(t1,[patterns,floor(101/disminuirbands)])')';
H1=reshape(mean(H(:,:,:,1),4),[M/disminuir*M/disminuir,floor(101/disminuirbands)]);%
H1 = normc(H1')';
distances = zeros(M/disminuir*M/disminuir,patterns);
for i=1:patterns
    t = t1(i,:)'; % compute the distance of the selected patterns with the whole matix
    distances(:,i) = H1*t;
end

[df,tf] = max(distances');
% tf stores the distribution of the repeated filters
Hs = cell(1,shots);
for k=1:shots
    Hs{k} = ((reshape(H(:,:,:,k),[M/disminuir*N/disminuir,floor(101/disminuirbands)]))')';
end
P = cell(1,patterns);
for i=1:patterns
    P{i}=zeros(shots,floor(101/disminuirbands));
    for j=1:shots
        P{i}(j,:) = mean(Hs{j}(tf==i,:),1); % extract the partitions 
    end
end
save('ready/H','H')
save('ready/tf','tf')
save('ready/P','P')
params.M = M;
params.N = N;
params.disminuir = disminuir;
params.patterns = patterns;
params.shots = shots;
save('ready/params','params')

