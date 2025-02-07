clc;clear;close all;

addpath(genpath(pwd));

I1 = imread('00000656.jpg');
I2 = imread('00000745.jpg');

tx = 15; ty=50;

if size(I1,3)>1
    img1 = I1;
else
    img1(:,:,1) = I1;img1(:,:,2) = I1; img1(:,:,3) = I1;
end

if size(I2,3)>1
    img2 = I2;
else
    img2(:,:,1) = I2;img2(:,:,2) = I2; img2(:,:,3) = I2;
end

%% feature extraction and display
[tar_feat,tar_desc, ref_feat, ref_desc] = sift_process(img1,img2);

if size(tar_feat,2)>4000
    indd = randperm(size(tar_feat,2),4000);
    tar_feat = tar_feat(:,indd);
    tar_desc = tar_desc(:,indd);
end

if size(ref_feat,2)>4000
    indd = randperm(size(ref_feat,2),4000);
    ref_feat = ref_feat(:,indd);
    ref_desc = ref_desc(:,indd);
end

height = max(size(img1,1),size(img2,1));
width = size(img1,2)+size(img2,2);

newfig=zeros(height, width,size(img1,3));
newfig(1:size(img1,1),1:size(img1,2),:) = img1;
newfig(1:size(img2,1),size(img1,2)+1:end,:)=img2;
newfig=uint8(newfig);
figure;
image(newfig);
axis image;axis off;

h_num = 20;
x1 = [0,width]';
y1 =  repmat(0:height/h_num:height,2,1);
line(x1,y1,'Color',[0.3,0.3,0.3]);
w_num = 20;
y2 = [0,height]';
x2 =  repmat(0:width/w_num:width,2,1);
line(x2,y2,'Color',[0.3,0.3,0.3]);
title('Griding');

figure;
width2 = max(size(img1,2),size(img2,2));
height2 = size(img1,1)+size(img2,1);

newfig2=zeros(height2, width2,size(img1,3));
newfig2(1:size(img1,1),1:size(img1,2),:) = img1;
newfig2(size(img1,1)+1:end,1:size(img2,2),:)=img2;
newfig2=uint8(newfig2);
image(newfig2); 
axis image; axis off;

if size(tar_feat,2)/2>2000
    pind_tf = randperm(size(tar_feat,2),2000);
else
    pind_tf = randperm(size(tar_feat,2),ceil(size(tar_feat,2)/2));
end
if size(ref_feat,2)/2>2000
    pind_rf = randperm(size(ref_feat,2),2000);
else
    pind_rf = randperm(size(ref_feat,2),ceil(size(ref_feat,2)/2));
end

f1 = tar_feat(:,pind_tf); f2 = ref_feat(:,pind_rf);
f2(2,:) = f2(2,:)+size(img1,1);
h1 = vl_plotframe(f1);
h2 = vl_plotframe(f2);
set(h1,'color','y','linewidth',0.5);
set(h2,'color','y','linewidth',0.5);
hold off
title('Feature Descriptors');

% setting of plot parameters
t_xd = 20; t_yd = 90; fontsize = 16;

%% feature matching --- my method GSLC

tic;
[~,seeds,candidates, matches] = GSLC(tar_feat, tar_desc, ref_feat, ref_desc, size(img1), size(img2));
toc;

% % seeds

figure;
image(newfig); 
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);hold on;

h_num = 20;
x1 = [0,width]';
y1 =  repmat(0:height/h_num:height,2,1);
line(x1,y1,'Color',[0.3,0.3,0.3]);
w_num = 20;
y2 = [0,height]';
x2 =  repmat(0:width/w_num:width,2,1);
line(x2,y2,'Color',[0.3,0.3,0.3]); hold on;

plot_matches = seeds';
plot(f1(1,plot_matches(1,:)), f1(2,plot_matches(1,:)), 'y.', f2(1,plot_matches(2,:)), f2(2,plot_matches(2,:)), 'y.', 'MarkerSize', 6);
hold on;
line([f1(1,plot_matches(1,:));f2(1,plot_matches(2,:))],[f1(2,plot_matches(1,:));f2(2,plot_matches(2,:))],'linewidth',0.1,'color','g');
hold on;
title('Seed Correspondences');

% % candidates

figure;
image(newfig); 
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);hold on;

h_num = 20;
x1 = [0,width]';
y1 =  repmat(0:height/h_num:height,2,1);
line(x1,y1,'Color',[0.3,0.3,0.3]);
w_num = 20;
y2 = [0,height]';
x2 =  repmat(0:width/w_num:width,2,1);
line(x2,y2,'Color',[0.3,0.3,0.3]);

if size(candidates,1)/3>800
    pind = randperm(size(candidates,1),800);
else
    pind = randperm(size(candidates,1),ceil(size(candidates,1)/3));
end

plot_matches = candidates(pind,:)';
plot(f1(1,plot_matches(1,:)), f1(2,plot_matches(1,:)), 'y.', f2(1,plot_matches(2,:)), f2(2,plot_matches(2,:)), 'y.', 'MarkerSize', 6.0);
hold on;
line([f1(1,plot_matches(1,:));f2(1,plot_matches(2,:))],[f1(2,plot_matches(1,:));f2(2,plot_matches(2,:))],'linewidth',0.1,'color','g');
hold on;
title('Candidate Correspondences');


figure;
image(newfig); 
axis image; axis off;
f1 = tar_feat; f2 = ref_feat;
f2(1,:) = ref_feat(1,:)+size(img1,2);hold on;

if size(matches,1)/3>800
    pind = randperm(size(matches,1),800);
else
    pind = randperm(size(matches,1),ceil(size(matches,1)/3));
end

plot_matches = matches(pind,:)';
plot(f1(1,plot_matches(1,:)), f1(2,plot_matches(1,:)), 'y.', f2(1,plot_matches(2,:)), f2(2,plot_matches(2,:)), 'y.', 'MarkerSize', 6.0);
hold on;
line([f1(1,plot_matches(1,:));f2(1,plot_matches(2,:))],[f1(2,plot_matches(1,:));f2(2,plot_matches(2,:))],'linewidth',1,'color','b');
hold on;
title('Matching Results');