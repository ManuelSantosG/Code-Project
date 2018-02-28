difference=[];
for i=1:81;
    difference=[difference abs(vv(:,index)'*vv(:,i))];
end
plot(difference)