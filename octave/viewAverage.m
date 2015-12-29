function []=viewAverage(dataset="Trace", nclass=4, lpath="../DATASETS")
close all;

file=strcat(lpath,'/',dataset,'.rel')
dsTr=load(file);
dsTr=relabelDS(dsTr);

[nTr,lTr]=size(dsTr);

for nc=1:nclass
  figure(nc)
  avr=zeros(1,lTr);
  deb=1;
  hold on;
  file=strcat('iDBA_centroid_', dataset, '_',num2str(nc), '.dat');
  cent1=load(file);
  plot(cent1(deb:lTr-1), 'linewidth', 2, 'color', 'b');
  file=strcat('iKDBA_centroid_', dataset, '_',num2str(nc), '.dat');  
  cent2=load(file);
  plot(cent2(deb:lTr-1), 'linewidth', 2, 'color', 'c');
  file=strcat('pKDTW_PWA_centroid_', dataset, '_',num2str(nc), '.dat');
  cent3=load(file);
  plot(cent3(deb:lTr-1), 'linewidth', 2, 'color', 'k');

  ii=0;
  for i=1:nTr
	if(dsTr(i,1)==nc || (nc==2 && dsTr(i,1)==-1))
		plot(dsTr(i,deb+1:lTr), 'r');
	end
  end
  for i=1:nTr
	if(dsTr(i,1)==nc || (nc==2 && dsTr(i,1)==-1))
		ii=ii+1;
		avr=avr+dsTr(i,deb:lTr);
	end
  end

plot(cent1(deb:lTr-1), 'linewidth', 2);
plot(cent2(deb:lTr-1), 'linewidth', 2, 'color', 'c');
plot(cent3(deb:lTr-1), 'linewidth', 2, 'color', 'k');
%plot(avr(deb+1:lTr)/ii,'-+k', 'linewidth', 2);
title(strcat('Class ',num2str(nc)));
legend('iDBA', 'iKDBA', 'pKDTW-PWA');
hold off;
end
end


function []=tt()
mat = load('mat-u1u2.dat');
[r,c]=size(mat);
matt=mat';
mat1=zeros(r,c);
mat2=zeros(r,c);
sc=sum(mat);
sr=sum(matt);
for i=1:r
  for j=1:c
     mat1(i,j)=mat(i,j)/sc(j);
     mat2(i,j)=mat(i,j)/sr(i);
  end
end
figure(1);
imagesc(mat1);
figure(2);
imagesc(mat2);
figure(3);
imagesc((mat1+mat2));
end

function [ds]=relabelDS(ds)
[n,l]=size(ds);
mn=min(ds(:,1));
mx=max(ds(:,1));

if(mn==-1 && mx==1)
  for i=1:n
  	if(ds(i,1)==-1)
 		ds(i,1)=1;
  	else
		ds(i,1)=2;
  	end
  end
else if (mn< 1 && mx>=0)
  for i=1:n
 	ds(i,1)=ds(i,1)+1; 
  end
end
end
end

