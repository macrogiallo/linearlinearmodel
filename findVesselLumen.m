function [M1]=findVesselLumen(I,M,num)
% manuscript: Deep Learning for quantification of basilar artery morphology
% using intracranial vessel wall MRI: a feasibility study

% input:
%   I: cross sectional vessel image
%   M: vessel outer wall boundary mask
%   num: number of probing directions
% ouput:
%   M1: vessel lumen (inner wall) mask


% num=12;% 12 directions
[s1,s2]=size(I);
[row,col,ix]=polarCoord_discrete(num,s1,s2);
px=[];
py=[];
idx=0;
M(M>0)=1;
for j=1:num
    y=I(ix(:,j));
    z=M(ix(:,j));
    ii=sum(z);

    if ii>=5
        [BAT,MD,TTP]=BAT_LLM_ridge(1:ii,double(y(1:ii)));
        i=MD;
    else
        [mx,imx]=max(y(1:ii));
        [mn,imn]=min(y(1:imx));
        mid=(mx+mn)/2;
        for i=imn:imx
            if y(i)>=mid
                break;
            end
        end
    end
    if i>0
        idx=idx+1;
        py(idx)=row(i,j);
        px(idx)=col(i,j);
    end
end
% disp(length(px))
M1=uint8(poly2mask(px,py,s1,s2));

function [row,col,ix] = polarCoord_discrete(n,s1,s2)
% create (discrete) polar coordinates indexing
% for fast radial sampling
% in n angular directions

m=floor(s1/2); % radius, assuming s1==s2
i0=m; % center of image
j0=m; % [s1/2,s1/2]
row=zeros(m,n); 
col=zeros(m,n);
dtheta=2*pi/n;
for j=1:n
    c=cos((j-1)*dtheta);
    s=sin((j-1)*dtheta);
    for i=1:m
        row(i,j)=i0+c*(i-1);
        col(i,j)=j0+s*(i-1);
    end
end
row=round(row); % discretized
col=round(col); % discretized
ix=sub2ind([s1,s2],row,col);

function [BAT,MD,TTP]=BAT_LLM_ridge(t,y)

% input:
%   t: time
%   y: contrast intensity
% ouput:
%   BAT: balus arrival time
%   MD: mid intensity location
%   TTP: time to peak

% LLM algorithm was adapted from:
% "An automatic approach for estimating bolus arrival time in dynamic 
% contrast MRI using piecewise continuous regression models"
% by LH Cheong, TS Koh, and Z Hou

% x2=length(y); % manual outer border
[~,x1]=max(y); % ymax at x1
TTP=t(x1); % converted to real time as t might be in a different time scale
[c1,c2]=size(y);
if c2>c1 % check if y is sized Nx1
    y=y';
end
C=y(1:x1);
x0=1;
err=realmax;% a big number
for i=2:x1-1 % linear search for the minimal fitting error
    X=ones(x1,2);
    X(:,2)=0;
    for j=i+1:x1
        X(j,2)=t(j)-t(i);
    end
%     B=inv(transpose(X)*X)*(transpose(X)*C); % find inv() explicitly
    B=(transpose(X)*X)\(transpose(X)*C); % '\' using gauss elimination
    D=C-X*B;
    e=norm(D);
    if e<err
        x0=i;
        err=e;
    end
end
BAT=t(x0); % refer to Scientific Reports 11(1) (2021) 479
md=(min(y(1:BAT))+y(TTP))/2; % mid intensity
for i=x0:x1
    if y(i)>=md
        MD=i;
        break;
    end
end
MD=t(MD); % converted to real time as t might be in a different time scale
