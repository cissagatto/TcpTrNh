##################################################################################################
# Test the Best Partition with CLUS                                                              #
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This code is free software: you can redistribute it and/or modify it under the terms of the    #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version. This code is distributed in the hope       #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################


FolderRoot = "~/TCP-TR-NH-Clus"
FolderScripts = "~/TCP-TR-NH-Clus/R"


##################################################################################################
# Runs for all datasets listed in the "datasets.csv" file                                        #
# n_dataset: number of the dataset in the "datasets.csv"                                         #
# number_cores: number of cores to paralell                                                      #
# number_folds: number of folds for cross validation                                             #
# delete: if you want, or not, to delete all folders and files generated                         #
##################################################################################################
execute <- function(ds, dataset_name, number_dataset,
                    number_folds, number_cores,
                    folderResults, diretorios){

  diretorios = diretorios


  FolderRoot = "~/TCP-TR-NH-Clus"
  FolderScripts = "~/TCP-TR-NH-Clus/R"

  if(number_cores == 0){
    cat("\n\n################################################################################################")
    cat("\n#Zero is a disallowed value for number_cores. Please choose a value greater than or equal to 1.#")
    cat("\n##################################################################################################\n\n")
  } else {
    cl <- parallel::makeCluster(number_cores)
    doParallel::registerDoParallel(cl)
    print(cl)

    if(number_cores==1){
      cat("\n\n###############################################################")
      cat("\n# ====> RUN: Running Sequentially!                              #")
      cat("\n#################################################################\n\n")
    } else {
      cat("\n\n###############################################################")
      cat("\n# ====> RUN: Running in parallel with ", number_cores, " cores! #")
      cat("\n#################################################################\n\n")
    }
  }

  retorno = list()

  setwd(FolderScripts)
  source("testClus.R")

  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Get dataset information: ", number_dataset, "          #")
  cat("\n#####################################################################\n\n")
  setwd(FolderRoot)
  datasets <- data.frame(read.csv("datasets-original.csv"))
  ds = datasets[number_dataset,]
  info = infoDataSet(ds)
  dataset_name = toString(ds$Name)


  cat("\n\n##################################################################")
  cat("\n# ====> RUN: VERIFYING THRESHOLDS                                  #")
  cat("\n#####################################################################\n\n")
  timeVerify = system.time(resTHP <- verifyingTresholds(ds,
                                                       dataset_name,
                                                       number_dataset,
                                                       number_folds,
                                                       number_cores,
                                                       folderResults,
                                                       diretorios))

  valid_tr = as.numeric(resTHP$valid_length)
  escolhidos = resTHP$tr_valid


  cat("\n\n###################################################################")
  cat("\n# ====> RUN: BUILD AND TEST PARTITIONS                              #")
  cat("\n#####################################################################\n\n")
  timeBuild = system.time(res <- buildAndTest(ds,
                                                 dataset_name,
                                                 number_dataset,
                                                 number_folds,
                                                 number_cores,
                                                 folderResults,
                                                 diretorios,
                                                 valid_tr))


  cat("\n\n###################################################################")
  cat("\n# ====> RU: Matrix Confusion                                        #")
  cat("\n#####################################################################\n\n")
  timeSplit = system.time(resGather <- juntaResultados(ds,
                                                       dataset_name,
                                                       number_dataset,
                                                       number_folds,
                                                       number_cores,
                                                       folderResults,
                                                       diretorios,
                                                       valid_tr))



  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Evaluation Fold                                        #")
  cat("\n#####################################################################\n\n")
  timeAvalia = system.time(resEval <- avaliaTest(ds,
                                                 dataset_name,
                                                 number_dataset,
                                                 number_folds,
                                                 number_cores,
                                                 folderResults,
                                                 diretorios,
                                                 valid_tr))



  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Gather Evaluation                                      #")
  cat("\n#####################################################################\n\n")
  timeGather = system.time(resGE <- juntaAvaliacoes(ds,
                                                    dataset_name,
                                                    number_dataset,
                                                    number_folds,
                                                    number_cores,
                                                    folderResults,
                                                    diretorios,
                                                    valid_tr))



  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Organize Evaluation                                    #")
  cat("\n#####################################################################\n\n")
  timeOrg = system.time(resGE <- organizaAvaliacoes(ds,
                                                    dataset_name,
                                                    number_dataset,
                                                    number_folds,
                                                    number_cores,
                                                    folderResults,
                                                    diretorios,
                                                    valid_tr))


  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Save Runtime                                           #")
  cat("\n#####################################################################\n\n")
  Runtime = rbind(timeVerify, timeBuild, timeBuild,
                  timeSplit, timeAvalia, timeGather,
                  timeGather,timeOrg)
  setwd(diretorios$folderTested)
  write.csv(Runtime, paste(dataset_name,
                           "-testPartition-Runtime.csv", sep=""),
            row.names = FALSE)

  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Stop Parallel                                          #")
  cat("\n#####################################################################\n\n")
  on.exit(stopCluster(cl))

  cat("\n\n###################################################################")
  cat("\n# ====> RUN: END                                                    #")
  cat("\n#####################################################################\n\n")
}

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
