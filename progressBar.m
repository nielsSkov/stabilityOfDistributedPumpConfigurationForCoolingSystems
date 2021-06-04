function progressBar( totalNrOfItrs, progressText )

persistent progressBarItrNr;
persistent progressBar_h;

if isempty(progressBar_h) || ~isgraphics(progressBar_h)
  progressBarItrNr = 1;
else
  progressBarItrNr = progressBarItrNr+1;
end

progress = sprintf( '%i/%i', progressBarItrNr,totalNrOfItrs);

if progressBarItrNr == 1
  progressBar_h = waitbar( progressBarItrNr/totalNrOfItrs, ...
                          [progressText, ' ', progress]         );
else
  waitbar( progressBarItrNr/totalNrOfItrs, progressBar_h, ...
          [progressText, progress]                        );
end

if progressBarItrNr == totalNrOfItrs
  pause(.1)
  close(progressBar_h)
end

pause(.001)

end

