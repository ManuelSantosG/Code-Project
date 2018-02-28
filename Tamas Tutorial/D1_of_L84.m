% 30.10.17. Demonstration of computing the information dimension of a
% fractal set living in the plane, e.g., the trajectory of a discrete time
% dynamical system, i.e., a map.

clear
close all

load XE

xsl = XE(:,1);
ysl = XE(:,2);

% Visualise the measure supported by the attractor
% Histogram bins
nb = 300; % number of bins in both directions, x and y
x_max = 3;
y_max = 3;
z_max = 1500;
edgesx = linspace(0,x_max,nb);
edgesy = linspace(0,y_max,nb);
% Color for histograms black, gary, gray. 
hcol = [1 1 1]*204/255; 

fh = figure;
set(fh,'Position',[100 100 600 500])

figure(fh)
hc = hist3([xsl, ysl],{edgesx, edgesy});
% Include in the same array the 1d histograms of the marginal distribution
hc2 = 0*hc;
hc2(1,:) = sum(hc,1);
hc2(1:end-1,end) = sum(hc(2:end,:),2);
% Truncate histogram bin counts
hc(hc>z_max) = z_max;
hc2(hc2>z_max) = z_max;
hand = bar3(5*hc,1);
    % 5* is an amplifier for visualisation; note that it makes no sense
    % anyway to use the same scale for the box counts of 2D and 1D
    % histograms, and so:
set(gca,'ZTickLabel',[]);
hold on
hand2 = bar3(hc2,1);
hold off
axis tight
xlabel('y'); ylabel('x'); zlabel('histogram box count') 
    % NOTICE the exchange of axis labels!!!

% 'bar3' has no such sensible axis-labeling as 'hist3'
% Solution from:
% http://stackoverflow.com/questions/17896848/plot-3d-histogram-using-bar3
dtick = round(nb/6); 
set(gca, 'xtick', 0.5:dtick:nb+0.5);
set(gca, 'ytick', 0.5:dtick:nb+0.5);
vec_bin_labels_x = ((1:dtick:nb+1)-1)*x_max/nb;
vec_bin_labels_y = ((1:dtick:nb+1)-1)*y_max/nb;
vec_string_bin_labels_x = ...
    reshape(cellstr(num2str(vec_bin_labels_x(:))), ...
    size(vec_bin_labels_x));
vec_string_bin_labels_y = ...
    reshape(cellstr(num2str(vec_bin_labels_y(:))), ...
    size(vec_bin_labels_y));
set(gca, 'xticklabel', vec_string_bin_labels_x);
set(gca, 'yticklabel', vec_string_bin_labels_y);

% We don't want to see bars with 0 height colored (corresponding to 0
% histogram bin counts)
% Solution from: 
% http://stackoverflow.com/questions/2050367/how-to-hide-zero-values-in-bar3-plot-in-matlab
for i2 = 1:numel(hand)
    index = logical(kron(hc(:,i2) == 0,ones(6,1)));
    zData = get(hand(i2),'ZData');
    zData(index,:) = nan;
    set(hand(i2),'ZData',zData,'EdgeColor','none','FaceColor','k');
end

for i2 = 1:numel(hand2)
  index = logical(kron(hc2(:,i2) == 0,ones(6,1)));
  zData = get(hand2(i2),'ZData');
  zData(index,:) = nan;
  set(hand2(i2),'ZData',zData,'FaceColor',hcol,'EdgeColor','none');
end


% Compute the information dimension of the measure
eps = exp(-linspace(1,7,100)); % resolution/size of box
I = 0*eps; % Number of boxes to fully cover attractor

for i = 1:length(eps)
    % Which 'x and y slots' is each point in?
    xb = floor(xsl/eps(i));
    yb = floor(ysl/eps(i));

    nx = max(xb) - min(xb) + 1; % number of x slots that span the object

    % Identify a box with a single number using the concept of digits
    xyb = xb + nx*yb;
    
    % Number of points in boxes
    Nb = diff(1+[0; find(diff(sort(xyb))); length(xyb)]);
    % Probability for boxes
    P = Nb/sum(Nb);
    
    % Information entropy of distribution
    I(i) = -sum(P.*log(P));
end

ord = I;
absc = log(1./eps);

figure; hold on
plot(absc,ord,'o','Color','magenta')
xlabel('ln(1/\epsilon)')
ylabel('I(\epsilon)')

% Fitting straight line
rg = 20:60; % fit range 

[p,S] = polyfit(absc(rg),ord(rg),1);
D1 = p(1);
f = polyval(p,absc(rg));
RMSE = sqrt(sum((f - ord(rg)).^2)/length(f));
title(['D1 = ' num2str(D1) ', RMSE = ' num2str(RMSE)])
hold on;
plot(absc,polyval(p,absc),'Color','k')
plot(absc(rg),f,'Color','g') % range where the fit took place