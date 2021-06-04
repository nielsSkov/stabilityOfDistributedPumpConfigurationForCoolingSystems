function saveCroppedPdf( figHandle, fileName )
	
	%check matlab version
	matVer = ver('MATLAB'); matVer = matVer.Release;
	releaseYear = sscanf(matVer,'(R%i');
	
	if matVer=="(R2020b)" || releaseYear >= 2021
		exportgraphics( figHandle, fileName )
	else	
		%set fig units to cm (default is pixel, not compatible with papersize)
		set(figHandle,'Units','centimeters');

		%get fig position to measure size of figure (used for cropping of pdf)
		pos = get(figHandle,'Position');

		%set figure paper properties
		set( figHandle, 'PaperPositionMode', 'Auto',  ...
										'PaperUnits', 'centimeters',  ...
										'PaperSize', [pos(3), pos(4)] )

		%print cropped figure to pdf
		print( figHandle, fileName, '-dpdf', '-r0')
	end
end