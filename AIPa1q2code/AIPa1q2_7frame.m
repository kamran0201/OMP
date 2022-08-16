clear;
close all;
clc;

rng('default');
% addpath('MMread');
[video,audio]=mmread('../cars.avi',1:7);
F=zeros(120,240,7,"single");
F(:,:,1)=rgb2gray(video.frames(1).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,1)));
F(:,:,2)=rgb2gray(video.frames(2).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,2)));
F(:,:,3)=rgb2gray(video.frames(3).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,3)));
F(:,:,4)=rgb2gray(video.frames(4).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,4)));
F(:,:,5)=rgb2gray(video.frames(5).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,5)));
F(:,:,6)=rgb2gray(video.frames(6).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,6)));
F(:,:,7)=rgb2gray(video.frames(7).cdata(169:288,113:352,:));
% figure();
% imshow(uint8(F(:,:,7)));

C=randi([0,1],120,240,7,'single');
N=2*randn(120,240,'single');
E=sum(C.*F,3)+N;
figure();
imshow(uint8(E/7));

RI=zeros(120,240,7,'single');
count=zeros(120,240,7,'single');
psi=single(dctmtx(448));
% psi=kron(eye(7),kron(dctmtx(8),dctmtx(8)));
% psi=kron(eye(7),dctmtx(64));


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
phi=reshape(C_patch,[],7);
phi=[diag(phi(:,1)),diag(phi(:,2)),diag(phi(:,3)),diag(phi(:,4)),diag(phi(:,5)),diag(phi(:,6)),diag(phi(:,7))];

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

x=zeros(448,1,'single');
x(T(:))=theta(:);

RI(u+1:u+8,v+1:v+8,:)=RI(u+1:u+8,v+1:v+8,:)+reshape(psi*x,8,8,7);
count(u+1:u+8,v+1:v+8,:)=count(u+1:u+8,v+1:v+8,:)+1;
% RI(u*8+1:u*8+8,v*8+1:v*8+8,:)=RI(u*8+1:u*8+8,v*8+1:v*8+8,:)+reshape(psi*x,8,8,7);
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
figure();
imshow(RI(:,:,6));
figure();
imshow(RI(:,:,7));

% rmse1=(norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)/(norm(F(:,:,1))^2);
% rmse2=(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)/(norm(F(:,:,2))^2);
% rmse3=(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2)/(norm(F(:,:,3))^2);
% rmse4=(norm(F(:,:,4)-single(RI(:,:,4)),'fro')^2)/(norm(F(:,:,4))^2);
% rmse5=(norm(F(:,:,5)-single(RI(:,:,5)),'fro')^2)/(norm(F(:,:,5))^2);
% rmse6=(norm(F(:,:,6)-single(RI(:,:,6)),'fro')^2)/(norm(F(:,:,6))^2);
% rmse7=(norm(F(:,:,7)-single(RI(:,:,7)),'fro')^2)/(norm(F(:,:,7))^2);
rmsetotal=((norm(F(:,:,1)-single(RI(:,:,1)),'fro')^2)+(norm(F(:,:,2)-single(RI(:,:,2)),'fro')^2)+(norm(F(:,:,3)-single(RI(:,:,3)),'fro')^2)+(norm(F(:,:,4)-single(RI(:,:,4)),'fro')^2)+(norm(F(:,:,5)-single(RI(:,:,5)),'fro')^2)+(norm(F(:,:,6)-single(RI(:,:,6)),'fro')^2)+(norm(F(:,:,7)-single(RI(:,:,7)),'fro')^2))/((norm(F(:,:,1))^2)+(norm(F(:,:,2))^2)+(norm(F(:,:,3))^2)+(norm(F(:,:,4))^2)+(norm(F(:,:,5))^2)+(norm(F(:,:,6))^2)+(norm(F(:,:,7))^2));

% fprintf('The RMSE for 1st frame is %0.4f\n', rmse1);
% fprintf('The RMSE for 2nd frame is %0.4f\n', rmse2);
% fprintf('The RMSE for 3rd frame is %0.4f\n', rmse3);
% fprintf('The RMSE for 4th frame is %0.4f\n', rmse4);
% fprintf('The RMSE for 5th frame is %0.4f\n', rmse5);
% fprintf('The RMSE for 6th frame is %0.4f\n', rmse6);
% fprintf('The RMSE for 7th frame is %0.4f\n', rmse7);
fprintf('The total RMSE is %0.4f\n', rmsetotal);





