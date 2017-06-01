############################################################
############ Simple cross-validation TEST 1
############################################################

############ Packages and functions

include("includeAllFiles.jl")

############ Test File

include("Z_FinalcBoP.jl")

########################################################################
########################################################################

############ Code for results

wdName = dirname(Base.source_path())
dirName = string(experimentName,"|",now())
mkdir(string(wdName,"/Results/",dirName))

### Loop on datasets
resultsDict = Dict()
resultsDict["nExtFolds"] = nExtFolds
resultsDict["nIntFolds"] = nIntFolds
resultsDict["nClassif"] = nClassif
for dataIndex in dataIndexVec

  println(string("/----------- DATASET: ",dataNames[dataIndex]," -----------/"))

  data = load(string(wdName,"/Datasets/",dataNames[dataIndex],".jld"))
  A = full(data["A"])
  A = (A + A')/2
  y = data["y"]
  n = length(y)

  ## Creating Splits
  extFoldsLabels,intFoldsLabels = createSplits(y,nExtFolds,nIntFolds)
  dataSetResultsDict = Dict()
  dataSetResultsDict["extFoldsLabels"] = extFoldsLabels
  dataSetResultsDict["intFoldsLabels"] = intFoldsLabels

  for classifIndex in 1:nClassif

    ### Loop on Methods
    for methodsIndex in methodsIndexVec

      println(string("METHOD #",methodsIndex," TEST"))

      tic()

      # Method extract
      preFuncToTest = allPreMethods[methodsIndex]
      preParValuesArray = allPreParameterValuesArray[methodsIndex]
      funcToTest = allMethods[methodsIndex]
      parValuesArray = allParameterValuesArray[methodsIndex]

      # Pre Method Handeling
      Karray = []
      if preParValuesArray != false
        for perPar in preParValuesArray
          K = preFuncToTest(A,perPar...)
          push!(Karray,K)
        end
      else
        K = preFuncToTest(A)
        push!(Karray,K)
      end


      ### Loop on external folds
      foldResultsValue = Float64[]
      bestPreParamVect = Float64[]
      bestParamVect = Float64[]
      for trainFold = 1:nExtFolds
        indexTrain = (extFoldsLabels .== trainFold)
        indexTest = (extFoldsLabels .!= trainFold)


        if preParValuesArray == false && parValuesArray == false
          bestPreParamVect = false
          bestParamVect = false
          labPrediction = funcToTest(K,y,indexTrain,indexTest)
        else
          #Internal CrossVal
          bestParams = bestParamCrossVal(Karray,y,funcToTest,preParValuesArray,parValuesArray,extFoldsLabels,trainFold,intFoldsLabels)

          if bestParams["bestParam"] == false
            labPrediction = funcToTest(Karray[bestParams["bestPreIndex"]],y,indexTrain,indexTest)
          else
            labPrediction = funcToTest(Karray[bestParams["bestPreIndex"]],y,indexTrain,indexTest,bestParams["bestParam"]...)
          end
          append!(bestPreParamVect,preParValuesArray[bestParams["bestPreIndex"]])
          append!(bestParamVect,bestParams["bestParam"])
        end

        totalNumber = sum(indexTest)
        nbrRight = sum(labPrediction .== y[indexTest])
        meanFoldResult = nbrRight/totalNumber

        println(string("External classification results: ",nbrRight,"/",totalNumber," (",round(meanFoldResult*100,3),"%)"))
        append!(foldResultsValue,meanFoldResult)
      end
      methodMeanStd = [mean(foldResultsValue),std(foldResultsValue)]
      elapseTime = toq()
      println(string("METHOD #",methodsIndex,"RESULTS: ",methodMeanStd*100)," (",elapseTime,"s)")
      methodResults = Dict("foldResultsValue" => foldResultsValue, "bestPreParamVect" => bestPreParamVect,"bestParamVect" => bestParamVect,"preParValuesArray" => preParValuesArray,"parValuesArray" => parValuesArray, "elapseTime" => elapseTime)
      save(string(wdName,"/Results/",dirName,"/",dataNames[dataIndex],"-C",classifIndex,"-Method",methodsIndex,"-res.jld"),methodResults)
      dataSetResultsDict[string("C",classifIndex,"-Method",methodsIndex)] = methodResults
    end
    resultsDict[string("C",classifIndex,"-M",dataNames[dataIndex])]= dataSetResultsDict
  end
end


save(string(wdName,"/Results/",dirName,"/Allresults.jld"),resultsDict)
resultsToCSV(dirName,dirName,dataIndexVec,dataNames,length(methodsIndexVec),nExtFolds,nClassif)
############
