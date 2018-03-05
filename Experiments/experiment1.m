a=1.4;
b=0.3;
tic
%Initial condition
x0=[0,0];

%Number of Iterations
n=10000;

%Sequence of images
Seq=zeros(2,n);
Seq(:,1)=x0;

for i=2:n
    v=Seq(:,i-1);
    [h1,h2]=henon(v(1),v(2),a,b);
    Seq(:,i)=[h1,h2];
end


plot(Seq(1,:),Seq(2,:),'.','MarkerSize',4)
title('Hénon Map')
hold on



v=Seq(:,n);
[h1,h2]=henon(v(1),v(2),a,b);
TmSeq=[Seq(:,2:n),[h1,h2]'];


% Seq=zeros(2,649);
% Seq(1,:)=XE([1:649],1);
% Seq(2,:)=XE([1:649],2);
% 
% TmSeq=zeros(649);
% TmSeq(1,:)=XE([2:650],1);
% TmSeq(2,:)=XE([2:650],2);
% 
% plot(Seq(1,:),Seq(2,:),'.','MarkerSize',4)
% hold on
% 
% n=649;
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
% Coarse Graining Boxes

nx=10;
ny=10;

xx=linspace(-1.5,1.5,nx);
yy=linspace(-0.5,0.5,ny);

% xx=linspace(0,3,nx);
% yy=linspace(0,2.5,ny);

%%%%%%%%%%%%%%%%%%%%
for i=1:nx
    for j=1:ny
        plot(xx(i),yy(j),'+','color','b')
        hold on
    end
end
%%%%%%%%%%%%%%%%%%%%
com_n=(nx-1)*(ny-1);

B=zeros(nx-1,ny-1,4);

for i=2:nx
    for j=2:ny
        B(i-1,j-1,1)=xx(i);
        B(i-1,j-1,2)=yy(j);
        B(i-1,j-1,3)=xx(i-1);
        B(i-1,j-1,4)=yy(j-1);
    end
end

Counter=zeros(nx-1,ny-1);
P=zeros(com_n,10);
for k=1:n
    Point=Seq(:,k);
    for i=1:nx-1
        for j=1:ny-1
            if Point(1)<=B(i,j,1) && Point(1)>B(i,j,3) && Point(2)<=B(i,j,2) && Point(2)>B(i,j,4)
                    Counter(i,j)=Counter(i,j)+1;
                    P(j+(ny-1)*mod(i-1,ny-1),Counter(i,j))=k;
            end
        end
    end
end

Bn=zeros(com_n,4);
for i=1:nx-1
    for j=1:ny-1
        Bn(j+(ny-1)*mod(i-1,ny-1),:)=B(i,j,:);
    end
end

% Counter=zeros(nx-1,ny-1);
% for k=1:n;
%     Point=Seq(:,k);
%     for i=1:nx-1;
%         for j=1:ny-1;
%             if Point(1)<=B(i,j,1) & Point(1)>B(i,j,3) & Point(2)<=B(i,j,2) & Point(2)>B(i,j,4);
%                     Counter(i,j)=Counter(i,j)+1;    
%             end
%         end
%     end
% end


Cn=zeros(com_n,1);
for i=1:nx-1
    for j=1:ny-1
        CCC=j+(ny-1)*mod(i-1,ny-1);
        Cn(CCC)=Counter(i,j);
    end
end



numerator=zeros(com_n,com_n);
for i=1:com_n
    for j=1:com_n
        if Cn(i)~=0
            for k=1:Cn(i)
                m=P(i,k);
                h=TmSeq(:,m);
                if h(1)<=Bn(j,1) && h(1)>Bn(j,3) && h(2)<=Bn(j,2) && h(2)>Bn(j,4)
                    numerator(i,j)=numerator(i,j)+1;
                end
            end
        end
    end
end

CCn=[];
iupa=[];
for i=1:com_n
    if Cn(i)~=0
        iupa=[iupa i];
        CCn=[CCn Cn(i)];
    end
end

new_n=length(CCn);
numerator_n=zeros(new_n);


q=1;
for i=1:com_n
    if Cn(i)~=0
        numerator_n(q,:)=numerator(i,iupa);
        q=q+1;
    end
end

TOp=zeros(new_n);
for i=1:new_n
    for j=1:new_n
        TOp(i,j)=numerator_n(i,j)/CCn(i);
    end
end


TrM=zeros(com_n,com_n);
for i=1:com_n
    for j=1:com_n
        if Cn(i)~=0
            TrM(i,j)=numerator(i,j)/Cn(i);
        end
    end
end

[vv,ll,ww]=eig(TOp);
lambda=diag(ll);
index=find(real(lambda<1+0.0005 & real(lambda>1-0.0005)));
Dv=ww(:,index);

hold off

Dv=(1/sum(Dv))*Dv;

figure(2)
plot(Dv)
title('Invariant Density')
% for i=1:com_n
%     if abs(Dv(i))>0
%         p1=0.5*(Bn(i,1)+Bn(i,3));
%         p2=0.5*(Bn(i,2)+Bn(i,4));
%         plot(p1,p2,'s','MarkerSize',10,'MarkerFaceColor','r')
%         %plot(p1,p2,'s','MarkerSize',2,'color','r')
%     end
% end
toc



