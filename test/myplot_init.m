function [h,handles,xx,yy,zz,holdState] = myplot_init(domain,nplotpts)

holdState = ishold();
% nplotpts = 51;
h = instantiateSlice3GUI();
handles = guihandles(h);
[xx, yy, zz] = meshgrid(linspace(domain(1), domain(2), nplotpts), ...
                        linspace(domain(3), domain(4), nplotpts), ...
                        linspace(domain(5), domain(6), nplotpts));

end