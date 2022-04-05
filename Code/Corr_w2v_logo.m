%%
clear
clc

load('RSA_w2vReal.mat');
load('RSA_logoReal.mat');

Zval_w2v = zeros(51,length(RSA_w2vReal_parsg_fm{1,1}));
for sub = 1:51
        temp = RSA_w2vReal_parsg_fm{sub,1};
        Zval_w2v(sub,:) = 0.5*(log((1+temp)./(1-temp)));       
end

Zval_logo = zeros(51,length(RSA_logoReal_parsg_fm{1,1}));
for sub = 1:51
        temp = RSA_logoReal_parsg_fm{sub,1};
        Zval_logo(sub,:) = 0.5*(log((1+temp)./(1-temp)));       
end

load('gfMask.mat');
rmask = load_untouch_nii('rVOT_all.nii');
mask_ref = spm_select('List',pwd,'^rgrayTPM.*\.nii$');
mask_hdrs = spm_vol(mask_ref);
[maskdata,~] = spm_read_vols(mask_hdrs);
maskINDs = find((maskdata>0.2).*(logical(rmask.img)~=0).*(gfMask~=0));

k1 = load_untouch_nii('exp2_w2v_RW_VOT_L_neg56.img');
k2 = load_untouch_nii('exp2_w2v_RW_VOT_L_neg34.img');
k3 = load_untouch_nii('exp2_logo_RW_VOT_L_neg54.img');
k4 = load_untouch_nii('exp2_logo_RW_VOT_L_neg30.img');
k5 = load_untouch_nii('exp2_logo_RW_VOT_R_neg52.img');

tp1 = find(k1.img(maskINDs)==1);rsa_fg1 = zeros(51,18);
for temp = 1:length(tp1)
    fg1 = Zval_w2v(:,tp1(temp));
    rsa_fg1(:,temp) = fg1;
end
tp2 = find(k2.img(maskINDs)==1);rsa_fg2 = zeros(51,length(tp2));
for temp = 1:length(tp2)
    fg2 = Zval_w2v(:,tp2(temp));
    rsa_fg2(:,temp) = fg2;
end
tp3 = find(k3.img(maskINDs)==1);rsa_fg3 = zeros(51,length(tp3));
for temp = 1:length(tp3)
    fg3 = Zval_logo(:,tp3(temp));
    rsa_fg3(:,temp) = fg3;
end
tp4 = find(k4.img(maskINDs)==1);rsa_fg4 = zeros(51,length(tp4));
for temp = 1:length(tp4)
    fg4 = Zval_logo(:,tp4(temp));
    rsa_fg4(:,temp) = fg4;
end
tp5 = find(k5.img(maskINDs)==1);rsa_fg5 = zeros(51,length(tp5));
for temp = 1:length(tp5)
    fg5 = Zval_logo(:,tp5(temp));
    rsa_fg5(:,temp) = fg5;
end

rsa_fg = zeros(51,5);
rsa_fg(:,1) = mean(rsa_fg1,2,'omitnan');%Sem-midFG-L
rsa_fg(:,2) = mean(rsa_fg2,2,'omitnan');%Sem-antFG-L
rsa_fg(:,3) = mean(rsa_fg3,2,'omitnan');%Logo-midFG-L
rsa_fg(:,4) = mean(rsa_fg4,2,'omitnan');%Logo-antFG-L
rsa_fg(:,5) = mean(rsa_fg5,2,'omitnan');%Logo-midFG-R

[r,p] = corr(rsa_fg,'type','Spearman');
save rsa_fg rsa_fg
%%
clear
clc

load('rsa_fg.mat');

myc = [155 126 194;24 189 180;255 174 185]/255;

subplot(2,4,1);%sementic corr(y=-56,y=-34)
scatter(rsa_fg(:,1),rsa_fg(:,2),'g','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,2),rsa_fg(:,1));
x = -.08:0.001:.1;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','-');hold on;
% title('Correlation between Semantic Representations');
text(.07,.07,'r = 0.26 marginal*','HorizontalAlignment','center');
text(.01,-.075,'Sem in left Middle FG','HorizontalAlignment','center');
text(-.11,.013,'Sem in left Anterior FG','HorizontalAlignment','center','rotation',90);

subplot(2,4,2);%logo corr(left y=-54,left y=-30)
scatter(rsa_fg(:,3),rsa_fg(:,4),'r','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,4),rsa_fg(:,3));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','--');hold on;
text(.07,.07,'r = 0.207 NS','HorizontalAlignment','center');
text(.01,-.075,'Logo in left Middle FG','HorizontalAlignment','center');
text(-.11,.013,'Logo in left Anterior FG','HorizontalAlignment','center','rotation',90);

subplot(2,4,3);%logo corr(left y=-54,right y=-52)
scatter(rsa_fg(:,3),rsa_fg(:,5),'r','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,5),rsa_fg(:,3));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','-');hold on;
text(.07,.07,'r = 0.485 ***','HorizontalAlignment','center');
% title('Correlation between Logo-grapheme Representations');
text(.01,-.075,'Logo in left Middle FG','HorizontalAlignment','center');
text(-.11,.013,'Logo in right Middle FG','HorizontalAlignment','center','rotation',90);

subplot(2,4,4);%logo corr(left y=-30,right y=-52)
scatter(rsa_fg(:,4),rsa_fg(:,5),'r','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,5),rsa_fg(:,3));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','-');hold on;
text(.07,.07,'r = 0.325 *','HorizontalAlignment','center');
text(.01,-.075,'Logo in left Anterior FG','HorizontalAlignment','center');
text(-.11,.013,'Logo in right Middle FG','HorizontalAlignment','center','rotation',90);

subplot(2,4,5);%anterior logo & anterior sem corr(left y=-30,left y=-34)
scatter(rsa_fg(:,2),rsa_fg(:,4),'b','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,4),rsa_fg(:,2));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','--');hold on;
text(.07,.07,'r = 0.195 NS','HorizontalAlignment','center');
text(.01,-.075,'Sem in left Anterior FG','HorizontalAlignment','center');%x
text(-.11,.013,'Logo in left Anterior FG','HorizontalAlignment','center','rotation',90);%y

subplot(2,4,6);%mid logo & mid sem corr(left y=-56,left y=-54)
scatter(rsa_fg(:,1),rsa_fg(:,3),'b','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,3),rsa_fg(:,1));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','--');hold on;
text(.07,.07,'r = 0.077 NS','HorizontalAlignment','center');
text(.01,-.075,'Sem in left Middle FG','HorizontalAlignment','center');%x
text(-.11,.013,'Logo in left Middle FG','HorizontalAlignment','center','rotation',90);%y

subplot(2,4,7);%anterior logo & mid sem
scatter(rsa_fg(:,1),rsa_fg(:,4),'b','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,~,rl,~] = regress(rsa_fg(:,4),rsa_fg(:,1));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','-');hold on;
text(.07,.07,'r = 0.284 *','HorizontalAlignment','center');
text(.01,-.075,'Sem in left Middle FG','HorizontalAlignment','center');%x
text(-.11,.013,'Logo in left Anterior FG','HorizontalAlignment','center','rotation',90);%y

subplot(2,4,8);%mid logo & anterior sem corr
scatter(rsa_fg(:,2),rsa_fg(:,3),'b','o','filled');hold on;
set(gca,'ytick',[ ]);set(gca,'ylim',[-.065,.09]);
set(gca,'xlim',[-.09,.11]);
set(gca,'xtick',[ ]);
[bl,bint,rl,rint] = regress(rsa_fg(:,3),rsa_fg(:,2));
x = -.11:0.001:.11;
y = bl*x + mean(rl);
plot(x,y,'k','Linestyle','--');hold on;
text(.07,.07,'r = 0.124 NS','HorizontalAlignment','center');
text(.01,-.075,'Sem in left Middle FG','HorizontalAlignment','center');%x
text(-.11,.013,'Logo in left Middle FG','HorizontalAlignment','center','rotation',90);%y

% 15 inch * 3 inch, 300dpi, font size 12.5 pound, lineweight 1.2 pound