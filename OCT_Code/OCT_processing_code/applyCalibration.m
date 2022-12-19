function pdIMAQCalibrated = applyCalibration(pdIMAQ, pdK, pnIndex, bMATLABInterp)

[nLineLength, nNumberLines] = size(pdIMAQ);

pdIMAQCalibrated = zeros([nLineLength, nNumberLines]);

if (bMATLABInterp)
    
    pnArray = (1 : nLineLength)';
    
    for nLine = 1 : nNumberLines
        pdIMAQCalibrated(:, nLine) = interp1(pdK, pdIMAQ(:, nLine), pnArray, 'spline');
    end
    
    clear pnArray

else
    
    for nLine = 1 : nNumberLines
        for nPoint = 1 : nLineLength
            nIndex = pnIndex(nPoint);
            
            fy1 = pdIMAQ(nIndex + 1, nLine);
            fy2 = pdIMAQ(nIndex + 2, nLine);
            fy3 = pdIMAQ(nIndex + 3, nLine);
            fy4 = pdIMAQ(nIndex + 4, nLine);
            
            fx1 = pdK(nIndex + 1) - nPoint;
            fx2 = pdK(nIndex + 2) - nPoint;
            fx3 = pdK(nIndex + 3) - nPoint;
            fx4 = pdK(nIndex + 4) - nPoint;
            
            fx1_2 = fx1 * fx1;
            fx2_2 = fx2 * fx2;
            fx3_2 = fx3 * fx3;
            fx4_2 = fx4 * fx4;
            
            fx1_3 = fx1 * fx1_2;
            fx2_3 = fx2 * fx2_2;
            fx3_3 = fx3 * fx3_2;
            fx4_3 = fx4 * fx4_2;
            
            f4 = fx1_3 * ((fx2_2 * fx3 * fy4 + fx3_2 * fx4 * fy2 + fx4_2 * fx2 * fy3) - (fx2_2 * fx4 * fy3 + fx3_2 * fx2 * fy4 + fx4_2 * fx3 * fy2)) - fx2_3 * ((fx1_2 * fx3 * fy4 + fx3_2 * fx4 * fy1 + fx4_2 * fx1 * fy3) - (fx1_2 * fx4 * fy3 + fx3_2 * fx1 * fy4 + fx4_2 * fx3 * fy1)) + fx3_3 * ((fx1_2 * fx2 * fy4 + fx2_2 * fx4 * fy1 + fx4_2 * fx1 * fy2) - (fx1_2 * fx4 * fy2 + fx2_2 * fx1 * fy4 + fx4_2 * fx2 * fy1)) - fx4_3 * ((fx1_2 * fx2 * fy3 + fx3_2 * fx1 * fy2 + fx2_2 * fx3 * fy1) - (fx1_2 * fx3 * fy2 + fx2_2 * fx1 * fy3 + fx3_2 * fx2 * fy1));
            f0 = fx1_3 * ((fx2_2 * fx3       + fx3_2 * fx4       + fx4_2 * fx2      ) - (fx2_2 * fx4       + fx3_2 * fx2       + fx4_2 * fx3      )) - fx2_3 * ((fx1_2 * fx3       + fx3_2 * fx4       + fx4_2 * fx1      ) - (fx1_2 * fx4       + fx3_2 * fx1       + fx4_2 * fx3      )) + fx3_3 * ((fx1_2 * fx2       + fx2_2 * fx4       + fx4_2 * fx1      ) - (fx1_2 * fx4       + fx2_2 * fx1       + fx4_2 * fx2      )) - fx4_3 * ((fx1_2 * fx2       + fx3_2 * fx1       + fx2_2 * fx3      ) - (fx1_2 * fx3       + fx2_2 * fx1       + fx3_2 * fx2      ));
            
            pdIMAQCalibrated(nPoint, nLine) = f4 / f0;
        end
    end
end

end