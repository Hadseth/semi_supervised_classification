############################################################
############ best paramCross val
############################################################

############ Function

bestParamCrossVal = function(Karray,y,funcToTest,preParValuesArray,parValuesArray,extFoldsLabels,extFold,intFoldsLabels)

  paramResultsValue = Float64[]

  for preParIndex in 1:length(preParValuesArray)
    prePar = preParValuesArray[preParIndex]
    K = Karray[preParIndex]

    for parIndex in 1:length(parValuesArray)
      par = parValuesArray[parIndex]

      foldResultsValue = Float64[]
      totalNumber = 0
      for intFold in sort(unique(intFoldsLabels))
        indexTrainInt = ((extFoldsLabels.==extFold) .& (intFoldsLabels.!=intFold))
        indexTestInt = ((extFoldsLabels.==extFold) .& (intFoldsLabels.==intFold))

        if par != false
          labPrediction = funcToTest(K,y,indexTrainInt,indexTestInt,par...)
        else
          labPrediction = funcToTest(K,y,indexTrainInt,indexTestInt)
        end

        append!(foldResultsValue,mean(labPrediction .== y[indexTestInt]))
        totalNumber += sum(indexTestInt)
      end

      meanFoldResult = mean(foldResultsValue)
      nbrRight = totalNumber * meanFoldResult
      println(string("Internal classification results: ",round(nbrRight),"/",totalNumber," (", round(meanFoldResult*100,3),"%)","| prePar= ",prePar," par=",par))
      append!(paramResultsValue,meanFoldResult)
    end
  end

  indMaximum = indmax(paramResultsValue)
  nP = length(parValuesArray)
  bestPreIndex = Int64(floor((indMaximum - 1)/nP)) + 1
  indPar = (indMaximum - 1)%nP + 1

  bestParam = parValuesArray[indPar]
  println(string("Best parameters: Pre=",preParValuesArray[bestPreIndex],", Par=",bestParam))

  return(Dict("bestPreIndex" => bestPreIndex, "bestParam" => bestParam))
end
