function saveCroppedPdf( figHandle, fileName )
	
	%set figure units to cm (default is pixel, not compatible with papersize)
	set(figHandle,'Units','centimeters');
	
	%get figure position to measure size of figure (used for cropping of pdf)
	pos = get(figHandle,'Position');
	
	%set figure paper properties
	set( figHandle, 'PaperPositionMode', 'Auto',  ...
	                'PaperUnits', 'centimeters',  ...
	                'PaperSize', [pos(3), pos(4)] )
	
	%print cropped figure to pdf
	print( figHandle, fileName, '-dpdf', '-r0')
end

% OBS: with MATLAB version R2021a and later, this is a cleaner solution:
% exportgraphics( gcf, fileName )