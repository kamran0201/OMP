clear;
close all;
clc;

rng('default');
% addpath('MMread');
[video,audio]=mmread('../cars.avi',1:3);
F=zeros(120,240,3,'single');
F(:,:,1)=rgb2gray(video.frames(1).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,1)));
F(:,:,2)=rgb2gray(video.frames(2).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,2)));
F(:,:,3)=rgb2gray(video.frames(3).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,3)));

C=randi([0,1],120,240,3,'single');
N=2*randn(120,240,'single');
E=sum(C.*F,3)+N;
figure();
imshow(uint8(E/3));

RI=zeros(120,240,3,'single');
count=zeros(120,240,3,'single');
psi=single(dctmtx(192));
% psi=kron(eye(3),kron(dctmtx(8),dctmtx(8)));
% psi=kron(ones(3),dctmtx(64));


for u=0:112
% tic;
for v=0:232
% for u=0:14
% tic;
% for v=0:29

E_patch=E(u+1:u+8,v+1:v+8);
C_patch=C(u+1:u+8,v+1:v+8,:);
% E_patch=E(u*8+1:u*8+8,v*8+1:v*8+8);
% C_patch=C(u*8+1:u*8+8,v*8+1:v*8+8,:);

y=reshape(E_patch,[],1);
phi=reshape(C_patch,[],3);
phi=[diag(phi(:,1)),diag(phi(:,2)),diag(phi(:,3))];

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

x=zeros(192,1,'single');
x(T(:))=theta(:);

RI(u+1:u+8,v+1:v+8,:)=RI(u+1:u+8,v+1:v+8,:)+reshape(psi*x,8,8,3);
count(u+1:u+8,v+1:v+8,:)=count(u+1:u+8,v+1:v+8,:)+1;
% RI(u*8+1:u*8+8,v*8+1:v*8+8,:)=RI(u*8+1:u*8+8,v*8+1:v*8+8,:)+reshape(psi*x,8,8,3);
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

% rmse1=(norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)/(norm(F(:,:,1))^2);
% rmse2=(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)/(norm(F(:,:,2))^2);
% rmse3=(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2)/(norm(F(:,:,3))^2);
rmsetotal=((norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)+(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)+(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2))/((norm(F(:,:,1))^2)+(norm(F(:,:,2))^2)+(norm(F(:,:,3))^2));

% fprintf('The RMSE for 1st frame is %0.4f\n', rmse1);
% fprintf('The RMSE for 2nd frame is %0.4f\n', rmse2);
% fprintf('The RMSE for 3rd frame is %0.4f\n', rmse3);
fprintf('The total RMSE is %0.4f\n', rmsetotal);










