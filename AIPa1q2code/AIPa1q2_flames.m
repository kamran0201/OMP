clear;
close all;
clc;

rng('default');
% addpath('MMread');
[video,audio]=mmread('../flame.avi',1:5);
F=zeros(288,352,5,"single");
F(:,:,1)=rgb2gray(video.frames(1).cdata);
% figure();
% imshow(uint8(F(:,:,1)));
F(:,:,2)=rgb2gray(video.frames(2).cdata);
% figure();
% imshow(uint8(F(:,:,2)));
F(:,:,3)=rgb2gray(video.frames(3).cdata);
% figure();
% imshow(uint8(F(:,:,3)));
F(:,:,4)=rgb2gray(video.frames(4).cdata);
% figure();
% imshow(uint8(F(:,:,4)));
F(:,:,5)=rgb2gray(video.frames(5).cdata);
% figure();
% imshow(uint8(F(:,:,5)));

C=randi([0,1],288,352,5,'single');
N=2*randn(288,352,'single');
E=sum(C.*F,3)+N;
figure();
imshow(uint8(E/5));

RI=zeros(288,352,5,'single');
count=zeros(288,352,5,'single');
psi=single(dctmtx(320));
% psi=kron(eye(5),kron(dctmtx(8),dctmtx(8)));
% psi=kron(ones(5),dctmtx(64));


for u=0:280
% tic;
for v=0:344
% for u=0:35
% tic;
% for v=0:43

E_patch=E(u+1:u+8,v+1:v+8);
C_patch=C(u+1:u+8,v+1:v+8,:);
% E_patch=E(u*8+1:u*8+8,v*8+1:v*8+8);
% C_patch=C(u*8+1:u*8+8,v*8+1:v*8+8,:);

y=reshape(E_patch,[],1);
phi=reshape(C_patch,[],5);
phi=[diag(phi(:,1)),diag(phi(:,2)),diag(phi(:,3)),diag(phi(:,4)),diag(phi(:,5))];

A=phi*psi;
A_vecnorm=vecnorm(A);
A_normalized=A./A_vecnorm;

r=y;
i=0;
theta=[];
T=single([]);
while(norm(r)>48)
    [maxi,j]=max(abs(r'*A_normalized));
    T=[T;j];
    i=i+1;
    A_j=A(:,T(:));
%     theta=(A_j'*A_j)\(A_j'*y);
    theta=pinv(A_j)*y;
    r=y-A_j*theta;
end

x=zeros(320,1,'single');
x(T(:))=theta(:);

RI(u+1:u+8,v+1:v+8,:)=RI(u+1:u+8,v+1:v+8,:)+reshape(psi*x,8,8,5);
count(u+1:u+8,v+1:v+8,:)=count(u+1:u+8,v+1:v+8,:)+1;
% RI(u*8+1:u*8+8,v*8+1:v*8+8,:)=RI(u*8+1:u*8+8,v*8+1:v*8+8,:)+reshape(psi*x,8,8,5);
% count(u*8+1:u*8+8,v*8+1:v*8+8,:)=count(u*8+1:u*8+8,v*8+1:v*8+8,:)+1;

end
% toc;
end


RI=uint8(RI./count);
figure();
imshow(RI(:,:,1));
figure();
imshow(RI(:,:,2));
figure();
imshow(RI(:,:,3));
figure();
imshow(RI(:,:,4));
figure();
imshow(RI(:,:,5));

% rmse1=(norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)/(norm(F(:,:,1))^2);
% rmse2=(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)/(norm(F(:,:,2))^2);
% rmse3=(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2)/(norm(F(:,:,3))^2);
% rmse4=(norm(F(:,:,4)-single(RI(:,:,4)),'fro')^2)/(norm(F(:,:,4))^2);
% rmse5=(norm(F(:,:,5)-single(RI(:,:,5)),'fro')^2)/(norm(F(:,:,5))^2);
rmsetotal=((norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)+(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)+(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2)+(norm(F(:,:,4)-single(RI(:,:,4)),'fro')^2)+(norm(F(:,:,5)-single(RI(:,:,5)),'fro')^2))/((norm(F(:,:,1))^2)+(norm(F(:,:,2))^2)+(norm(F(:,:,3))^2)+(norm(F(:,:,4))^2)+(norm(F(:,:,5))^2));

% fprintf('The RMSE for 1st frame is %0.4f\n', rmse1);
% fprintf('The RMSE for 2nd frame is %0.4f\n', rmse2);
% fprintf('The RMSE for 3rd frame is %0.4f\n', rmse3);
% fprintf('The RMSE for 4th frame is %0.4f\n', rmse4);
% fprintf('The RMSE for 5th frame is %0.4f\n', rmse5);
fprintf('The total RMSE is %0.4f\n', rmsetotal);





