function [h,handles] = myplot(domain,nplotpts,h,handles,xx,yy,zz,holdState,v,nSteps)

if isreal(v)
    [row,col,tube] = ind2sub(size(v), find(v(:) == max(v(:)), 1, 'last'));
else
    [row, col, tube] = ind2sub(size(v), find(abs(v(:)) == max(abs(v(:))), 1, 'last'));    
end
xslice = xx(row,col,tube); 
yslice = yy(row,col,tube); 
zslice = zz(row,col,tube); 

set(handles.xSlider, 'Min', domain(1));
set(handles.xSlider, 'Max', domain(2));
set(handles.xSlider, 'Value', xslice);

set(handles.ySlider, 'Min', domain(3));
set(handles.ySlider, 'Max', domain(4));
set(handles.ySlider, 'Value', yslice);

set(handles.zSlider, 'Min', domain(5));
set(handles.zSlider, 'Max', domain(6));
set(handles.zSlider, 'Value', zslice);

% nSteps = 15; % number of slices allowed
set(handles.xSlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.ySlider, 'SliderStep', [1/nSteps , 1 ]);
set(handles.zSlider, 'SliderStep', [1/nSteps , 1 ]);

% Choose default command line output for the slice command:
handles.xx = xx;
handles.yy = yy;
handles.zz = zz;
handles.xslice = xslice;
handles.yslice = yslice;
handles.zslice = zslice;
handles.v = v;

if isreal(v)
    slice(xx, yy, zz, v, xslice, yslice, zslice)
    shading interp
    colorbar
else
    hh = slice(xx, yy, zz, angle(-v), xslice, yslice, zslice); 
    set(hh, 'EdgeColor','none')
    caxis([-pi pi]),
    colormap('hsv')
    axis('equal')     
end

hold on

hold off

axis equal
xlim(domain(1:2))
ylim(domain(3:4))
zlim(domain(5:6))

% Update handles structure
guidata(h, handles);
handles.output = handles.xSlider;

if ( ~holdState )
    hold off
end

% Force the figure to clear when another plot is drawn on it so that GUI
% widgets don't linger.  (NB:  This property needs to be reset to 'add' every
% time we change the plot using a slider; otherwise, the slider movement will
% itself clear the figure, which is not what we want.)
set(h, 'NextPlot', 'replacechildren');


end