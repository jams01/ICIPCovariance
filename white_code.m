random = 'True';
shots = 8;

M = 512;
N = 512;
bands=31;
disminuirbands=3.25;

start1=513;
start2=64;

folder_path = ['codigos/',random,'_',num2str(shots),'/'];

offset = 7;
offsety = 0;
H = zeros(M,N,101);

for i = 1:shots
    code_path = [folder_path,num2str(i),'_aperture',num2str(shots),'rand=',lower(random),'.mat'];
    load(code_path);
    h = X(start1+offsety:start1+offsety+M-1,start2-offset:start2-offset+N-1,:);
    H = H + h;
   
end

H = H;
Hr = zeros(M,N,bands);
index=1;
bandsall =450:2:650;
bandsdisc = zeros(1, bands);
load('disper')

for i=1:31
    Hr(:,:,index) = sum(H(:,:,disper==i),3)./sum(disper==i);
    [v,p]=find(disper==i);
    bandsdisc(i)=mean(bandsall(p));
     index=index+1;
end

% normalizar por el maimo banda a banda


