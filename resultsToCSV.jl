############################################################
############################################################
############################################################
############################################################

using JLD
using DataFrames
using Distributions
using HypothesisTests

############################################################
############################################################
############### Result to CSV

resultsToCSV = function(folderName,resName,dataIndex,dataNames,nMethods,n_samp,n_class)

  wdName = dirname(Base.source_path())
  if !isdir(string(wdName,"/ResultsCSV/",resName))
    mkdir(string(wdName,"/ResultsCSV/",resName))
  end
  dataF = DataFrame()
  dataF[:dataNames]= dataNames[dataIndex]

  methodNames = [string("Method",i) for i in 1:nMethods]

  allScore = Array{Float64}(length(dataIndex),0)
  for met in methodNames
    meanVec = Float64[]
    sdVec = Float64[]
    methScores = Float64[]
    for set in dataIndex
      totalFoldRes = Float64[]
      for classifInd in 1:n_class
        res = load(string(wdName,"/Results/",folderName,"/",dataNames[set],"-C",classifInd,"-",met,"-res.jld"))
        append!(totalFoldRes,res["foldResultsValue"])
        #append!(methScores,mean(res["foldResultsValue"]*100/90))
      end
      append!(meanVec,mean(totalFoldRes)*100)
      append!(methScores,mean(totalFoldRes))
      append!(sdVec,std(totalFoldRes)/sqrt(length(totalFoldRes))*100)
    end
    allScore = [allScore methScores]
    dataF[Symbol(string(met,"M"))] = meanVec
    dataF[Symbol(string(met,"SD"))] = sdVec
  end
  writetable(string(wdName,"/ResultsCSV/",resName,"/Res_.csv"), dataF,separator = ';')

  rankData = DataFrame()
  rankData[:methodNames]= methodNames
  for datInd in 1:length(dataIndex)
    rankData[Symbol(string("F",datInd))] = sortperm(sortperm(allScore[datInd,:]))
  end
  meanR = mean(Matrix(rankData[:,2:end]),2)
  stdR = std(Matrix(rankData[:,2:end]),2)/sqrt(size(rankData)[2]-1)
  rankData[:meanrank] = meanR[:]
  rankData[:stdmeanrank] = stdR[:]
  writetable(string(wdName,"/ResultsCSV/",resName,"/Res_MeanRank.csv"), rankData,separator = ';')

  normD = Normal()
  (nData,nMethods) = size(allScore)
  wilcoxonMatrix = zeros(Float64,nMethods,nMethods)
  for m1 in 1:nMethods
    for m2 in 1:nMethods
      m1Vec = allScore[:,m1]
      m2Vec = allScore[:,m2]
      pValue = pvalue(SignedRankTest(m1Vec, m2Vec),tail=:right)
      wilcoxonMatrix[m1,m2] = pValue
    end
  end
  wilcoxonMatrixOrdered = wilcoxonMatrix[sortperm(meanR[:],rev=true),sortperm(meanR[:],rev=true)]

  writetable(string(wdName,"/ResultsCSV/",resName,"/Res_wilcoxonMatrix.csv"), DataFrame(wilcoxonMatrixOrdered),separator = ';')

end
