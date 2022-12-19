function pdMask = calculateMask(nLineLength, nLeft, nRight, nRound)

pdMask = zeros([nLineLength 1]);


if (nRound < 1)
    nRound = 1;
end
if (nRound > (nLineLength / 4))
    nRound = nLineLength / 4;
end
if (nLeft < nRound)
    nLeft = nRound + 1;
end
if (nLeft > nLineLength - nRound)
    nLeft = nLineLength - nRound;
end
if (nRight < nRound + 1)
    nRight = nRound + 1;
end
if (nRight > nLineLength - nRound - 1)
    nRight = nLineLength - nRound - 1;
end
if (nLeft >= nRight)
    nRight = nLeft + 1;
end


for i = nLeft - nRound : nLeft
    pdMask(i) = 0.5 * (1 + cos(pi * (nLeft - i) / nRound));
end
for i = nLeft : nRight
    pdMask(i) = 1;
end
for i = nRight : nRight + nRound
    pdMask(i) = 0.5 * (1 + cos(pi * (i - nRight) / nRound));
end

end