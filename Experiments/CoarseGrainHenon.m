
res=0.1;
a=1.4;
b=0.3;

xx=[0:res:1];
yy=[0:res:1];

nx=length(xx);
ny=length(yy);

domain=zeros(2,nx*ny);
k=0;
for i=1:nx;
    for j=1:ny;
        k=j+(mod(i-1,nx))*ny
        domain(:,k)=[xx(i),yy(j)];
    end
end

for k=1:nx*ny;
    [domain(1,k),domain(2,k)]=henon(domain(1,k),domain(2,k),a,b);
end

for i=1:nx*ny;
    plot(domain(1,i),domain(2,i),'*','MarkerSize',4)
    hold on
end

% v=zeros(2,nx*ny);
% for i=1:nx;
%     for j=1:ny;
%         v(:,i)=henon(xx(i),yy(i),a,b);
%     end
% end

% for i=1:nx*ny;
%     plot(v(1,i),v(2,i),'*','MarkerSize',4)
%     hold on
% end

