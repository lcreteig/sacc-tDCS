function countDownTimer(windowPtr,xCoord,yCoord,startTime)

counter = startTime;
while counter > 0
    textDisp = sprintf('Task starts in %i', counter);
    DrawFormattedText(windowPtr, textDisp, xCoord, yCoord);
    Screen('Flip', windowPtr);
    WaitSecs(1);
    counter = counter-1;
end