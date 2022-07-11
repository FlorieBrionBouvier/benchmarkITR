# This script implement the functions used to compute the crossfited estimator
# and non crossfited estimator.
# User should use the function "compute_method"

#### Initialization ####

# setwd("my work directory")
rm(list = ls())
graphics.off()
set.seed(123)


# library import

library(rms) # cubic spline
library(magrittr) # pipeline
library(tidyverse) # data.frame manipulation
library(randomForest) # radom forest model


# import dataset

# ist <- read.csv("C:/Users/etien/Documents/cours/Benchmark ITR/ist_partitioned.csv")
# # ist <- read.csv("C:\\Users\\Florie BRION-BOUVIER\\Documents\\These\\Benchmark_ITR\\ist_partitioned.csv")
# df <- ist


#### Pre Processing ####

# # Changing the reference label of some categorical variables to match Collins
# df <- df %>% mutate(Conscious = relevel(factor(Conscious), "Fully alert"),
#                     StrokeType = relevel(factor(StrokeType), "PACS"),
#                     region = relevel(factor(region), "Europe"))
# 
# # train and test dataset
# train <- filter(df, group == "Train")
# test <- filter(df, group == "Test")


#### Functions ####

.idx_na_col = function(l) {
  # get the index of the columns containing NA values
  apply(l, 2,
        function(col) any(is.na(col))
        ) %>% unlist %>% which
}


add_dependencies = function(Balancing_Method) {
  # for each nuisance function of eache method, add the required number of split
  # and the nuisance functions needed to be called to compute the nuisance function
  
  lapply(
    Balancing_Method,
    function(met) {
      # for each nuisance function get the names of the other nuisance function
      # required to compute the nuisance function
      args_name = lapply(
        met[["nuisance_fct"]],
        function(nui){
          formalArgs( nui[["train"]] ) %>%
            `[`( . %in% names(met[["nuisance_fct"]]) ) %>%
            { if (length(.)==0) NULL else . }
          })
      
      # if a nuisance function do not require other nuisance function then it
      # only requires one split
      nbr_required_split = sapply(args_name, function(name) ifelse(is.null(name), 1, NA))
      
      # if a function nuisance do not require other nuisance function then it
      # does not have any dependencies
      dependencies = lapply(args_name, function(name) if (is.null(name)) list() else NA)
      
      
      while(any(is.na(nbr_required_split))) {
        # get the names of nuisance function  for which 'nbr_required_split' and
        # 'dependencies' are still not computed
        to_compute = nbr_required_split %>%
          is.na %>%
          nbr_required_split[.] %>%
          names
        
        for (name in to_compute) {
          # the number of split required by a nuisance function is the sum of
          # required number split of its direct dependencies + 1
          nbr_required_split[name] = args_name[[name]] %>%
            nbr_required_split[.] %>%
            sum(1)
          
          # the dependencies of a nuisance function is the union of the
          # dependencies's of the nuisance function it requires to be computed
          dependencies[[name]] = args_name[[name]] %>%
            met[["nuisance_fct"]][.] %>%
            lapply(`[[`, "train") %>%
            c(dependencies[[args_name[[name]]]]) %>%
            unique %>%
            { if (any(is.na(.))) NA else . }
        }
      }
      
      # save the values
      for (nui_name in names(met[["nuisance_fct"]])) {
        met[["nuisance_fct"]][[nui_name]]["nbr required split"] = nbr_required_split[nui_name]
        met[["nuisance_fct"]][[nui_name]][["dependencies"]] = dependencies[[nui_name]]
      }
      return(met)
    }
  )
}

add_unique_nuisance = function(Balancing_Method) {
  # create the skeleton for 'unqiue_nuisance', a list containing the minimum
  # number of nuisance function needed to be computed to crossfit the methods in
  # 'Balancing_Method'.
  #
  # for each nuisance function, 'unique_nuisance' contains :
  #   - the nuisance function;
  #   - the name given to this nuisance function by the 1st method in Balancing_Method
  #     using the nuisance function;
  #   - the names of the other nuisance functions required to compute the nuisance
  #     function (the names are also retrived from the 1st method calling the nuisance function),
  #   - the index in 'unique_nuisance' of the other nuisance functions required to compute
  #     the nuisance function;
  #   - a boolean equal to TRUE if the nuisance function should predict the values in
  #     the other splits;
  #   - a boolean equal to TRUE if the nuisance function should predict the values in test_dataset;
  #   - the indexes of the methods requiring the nuisance function.
  
  # check if the number of split for each method is the same
  sapply(Balancing_Method, `[[`, "nbr_split") %>%
    { any(.[1] != .) } %>%
    if (.) stop("'add_unique_nuisance' should only be used for 'Balancing_Method'
                containing method with equal 'nbr_split'")
  
  nbr_split = Balancing_Method[[1]][["nbr_split"]]
  
  # add dependencies if missing
  sapply(
    Balancing_Method,
    function (met) {
      sapply(met[["nuisance_fct"]], `[[`, "nbr required split") %>%
        sapply(is.null) %>% any
    }) %>%
    any %>%
    if (.) Balancing_Method <<- add_dependencies(Balancing_Method)
  
  # create unique_nuisance
  unique_nuisance = lapply(Balancing_Method, `[[`, "nuisance_fct") %>% 
    # concatenate the nuisance functions
    do.call(c,.) %>% 
    # remove duplicate
    unique %>% 
    # order the list
    `[`(sort.list(sapply(., `[[`, "nbr required split")))
  
  
  Balancing_Method[["unique_nuisance"]] = lapply(
    unique_nuisance,
    function(nui) {
        
      
      # get the name of 'nui' and a method calling 'nui'.
      # get the index of the method requiring 'nui'
      met_idx = sapply(
        seq_along(Balancing_Method),
        function(i) {
          match(list(nui), Balancing_Method[[i]][["nuisance_fct"]]) %>%
            sapply(is.na) %>% all(.)
        }) %>% `!` %>% which
      met = Balancing_Method[[first(met_idx)]]
      nui_name = match(list(nui), met[["nuisance_fct"]]) %>%
        names(met[["nuisance_fct"]])[.] %>% first
      
      
      # should 'nui' predict the value of the test dataset
      predict_test_dataset = sapply(
        Balancing_Method,
        function(met) {
          if (!(nui_name %in% formalArgs(met[["method"]]))) return(FALSE)
          return(list(nui) %in% met[["nuisance_fct"]])
        }) %>% any
      
      
      # should 'nui' predict on the other splits
      predict_split = sapply(
        unique_nuisance,
        function(other_nui) {
          if (!(nui_name %in% formalArgs(other_nui[["train"]]))) return(FALSE)
          return(nui["train"] %in% other_nui[["dependencies"]])
        }) %>% any
      
      
      # names of the required nuisance functions for 'nui' 
      names_required_nui = names( met[["nuisance_fct"]] ) %>%
        `[`( . %in% formalArgs(nui[["train"]]) ) %>%
        { if (length(.) == 0) NULL else . }
      
      # idx of the required nuisance function in unique_nuisance
      nui_unique_nuisance_idx = if (is.null(names_required_nui)) NULL else sapply(
        seq_along(names_required_nui),
        function(j) {
          # get the name of the required nuisance function
          met[["nuisance_fct"]][names_required_nui][[j]][c("train", "dependencies")] %>%
            list %>%
            match(lapply(unique_nuisance, `[`, c("train", "dependencies")))
          })
      
      
      lapply(
        1:nbr_split,
        function(i) {
          
          if ( length(nui[["dependencies"]]) == 0 ) return(NULL)
          
          # get the index of the split used to train the required nuisance
          # functions
          names_required_nui %>%
            { met[["nuisance_fct"]][.] } %>%
            sapply(`[[`, "nbr required split") %>%
            c(0,.) %>% cumsum %>% .[-length(.)] %>%
            { i - . - 1 } %>%
            { (.-1) %% nbr_split + 1} %>%
            list("nui_train_split_idx" = .)
        }) %>%
        c(nui,
          list("nui_name" = nui_name,
               "ignore" = FALSE,
               "met_idx" = met_idx,
               "predict_split" = predict_split,
               "predict_test_dataset" = predict_test_dataset,
               "names_required_nui" = names_required_nui,
               "nui_unique_nuisance_idx" = nui_unique_nuisance_idx)
          )
    })
  return(Balancing_Method)
}


add_method_info = function(Balancing_Method) {
  # creat 'method_info', a list containing useful information (such as the index
  # in 'unique_nuisance' of the nuisance function needed to compute a method) to
  # crossfit a method in 'Balancing_Method'.
  
  # check if 'unique_nuisance' is already in 'Balancing_Method'
  if (is.null(Balancing_Method[["unique_nuisance"]]))
    Balancing_Method %<>% add_unique_nuisance
  
  unique_nuisance = Balancing_Method[["unique_nuisance"]]
  
  nbr_split = Balancing_Method[[1]][["nbr_split"]]
  Balancing_Method[["method_info"]] = lapply(
    Balancing_Method[-length(Balancing_Method)],
    function(met) {
      
      # names of the nuisance functions needed to be passed on as an arg to
      # met[["method"]]
      names_required_nui = names(met[["nuisance_fct"]]) %>%
        `[`( . %in% formalArgs(met[["method"]]) ) %>%
        { if (length(.) == 0) NULL else . }
      
      # index in 'unique_nuisance' of the nuisance functions needed to be passed
      # on as an arg to met[["method"]]
      nui_unique_nuisance_idx = if (is.null(names_required_nui)) NULL else sapply(
        seq_along(names_required_nui), 
        function(i) {
          met[["nuisance_fct"]][names_required_nui][[i]][c("train", "dependencies")] %>%
            list %>%
            match( lapply(unique_nuisance,`[`, c("train", "dependencies")) )
        })
      
      
      lapply(
        1:nbr_split,
        function(split_idx) {
          # if met does
          if (is.null(names_required_nui)) return(NULL)
          
          # get the index of the split used to train the required nuisance
          # functions
          met[["nuisance_fct"]][names_required_nui] %>%
            sapply(`[[`, "nbr required split") %>%
            c(0, .) %>% cumsum %>% `[`( -length(.) ) %>%
            { split_idx - . - 1 } %>%
            { (.-1) %% nbr_split + 1} %>%
            list("nui_train_split_idx" = .)
          }) %>%
        c(list("names_required_nui" = names_required_nui,
               "nui_unique_nuisance_idx" = nui_unique_nuisance_idx)
          )
    })
  return(Balancing_Method)
}


subset_Balancing_Method = function (Balancing_Method, idx) {
  to_ignore = sapply(
    lapply(Balancing_Method[["unique_nuisance"]], `[[`, "met_idx"),
    function(i) i %in% idx %>% any %>% `!`
  ) %>% which
  
  for (i in to_ignore) Balancing_Method[["unique_nuisance"]][[i]][["ignore"]] = TRUE
  
  c(Balancing_Method[idx],
    list("method_info" = Balancing_Method[["method_info"]][idx],
         "unique_nuisance" = Balancing_Method[["unique_nuisance"]])
  )
}


simple_crossfit = function(X_tr, Y_tr, Z_tr, X_te, Y_te, Z_te,
                           Balancing_Method, nbr_iter, nbr_split) {
  n = nrow(X_tr) # number of observations in train dataset
  p = ncol(X_tr) # number of covariates
  
  if (is.null(Balancing_Method[["method_info"]]))
    Balancing_Method %<>% add_method_info
  
  unique_nuisance = Balancing_Method[["unique_nuisance"]]
  method_info = Balancing_Method[["method_info"]]
  Balancing_Method = Balancing_Method[-length(Balancing_Method) + 0:1]
  
  test_dataset = list(
    "X" = X_te,
    "Y" = Y_te,
    "Z" = Z_te)
  
  # split the dataset "nbr_iter" times
  
  replicate(
    n = nbr_iter,
    simplify = FALSE,
    expr =
    {
      # split the data evenly without overlap
      splited_data = cbind(X_tr, Y_tr, Z_tr) %>%
        # shuffle the rows
        .[sample(n), ] %>%
        # split [X_tr, Y_tr, Z_tr] into "nbr_split" splits
        split( rep_len(1:nbr_split, n) ) %>% 
        # unbind [X_tr, Y_tr, Z_tr]
        lapply( function(data) list(
          "X" = data[, 1:p],
          "Y" = data[, p+1],
          "Z" = data[, p+2])
        )
      
      
      # train the models on each split and predict the values on the other splits
      for (nui_idx in seq_along(unique_nuisance)) {
        
        nui = unique_nuisance[[nui_idx]]
        
        if (nui[["ignore"]]) next
        
        unique_nuisance[[nui_idx]] = lapply(
          1:nbr_split,
          function(i) {
            
            out = nui[i]
            
            # get the predicted values on the split i of the required nuisance
            # functions trained on the split nui_train_split_idx
            required_nui_input =  lapply(
              seq_along(nui[["names_required_nui"]]),
              function(j) {
                # get the prediction on the split i
                unique_nuisance %>%
                  `[[`(nui[["nui_unique_nuisance_idx"]][j]) %>%
                  `[[`(nui[[i]][["nui_train_split_idx"]][j]) %>%
                  `[[`("predicted.value") %>%
                  `[`(i)
              }) %>%
              `names<-`( nui[["names_required_nui"]] ) %>%
              do.call(c, .)
            
            
            # train the nuisance function on the split i
            out[["trained.model"]] = splited_data[[i]] %>%
              c(required_nui_input) %>%
              # remove unused args
              `[`( names(.) %in% formalArgs(nui[["train"]]) ) %>%
              { try( do.call(nui[["train"]], .) ) } %>%
              # return NA if "nui" fail
              { if (any(class(.) == "try-error")) NA else . }
            
            
            # predicted the values on the other splits and test_dataset
            out[["predicted.value"]] = lapply(
              1:(nbr_split+1),
              function(j) {
                # if there is no need to predict on the splits
                if (j <= nbr_split & !nui[["predict_split"]]) return(NULL)
                # if there is no need to predict on the test dataset
                if (j == nbr_split & !nui[["predict_test_dataset"]]) return(NULL)
                # do not predict on the split used to train the model
                if (j==i) return(NULL)
                
                # do not try to predict if the training failed
                if ( all(is.na(out[["trained.model"]])) ) return(NA)
                if ( !is.null(out[["trained.model"]][["fail"]]) )
                  if ( out[["trained.model"]][["fail"]] ) return(NA)
                
                # predict on the other splits and the testing data set
                (if (j == nbr_split+1) test_dataset else splited_data[[j]]) %>%
                  c(out["trained.model"]) %>%
                  # remove unused args
                  .[names(.) %in% formalArgs(nui[["predict"]])] %>%
                  # preidict the values in the current split
                  { try(do.call(nui[["predict"]], .)) } %>%
                  # return NA if "mnui$predict" fail
                  { if (any(class(.) == "try-error")) NA else . }
              }) %>%
              # rename the lists containing the predictions to 1, 2, ..., nbr_split, "test_dataset"
              `names<-`( c(rep("",nbr_split),"test_dataset") )
            
            return(out)
          }) %>%
          c(nui[-(1:nbr_split)])
      }
      
      
      # crossfit each method
      lapply(
        seq_along(Balancing_Method),
        function(met_idx) {
          met = Balancing_Method[[met_idx]]
          output = lapply(
            1:nbr_split,
            function(split_idx) {
              # determine the index of the splits used to train
              # the nuisance functions for the current method
              nui_split_idx = method_info[[met_idx]][[split_idx]][["nui_train_split_idx"]]
              
              
              # create a list containing the arguments needed to compute
              # met[["method"]]
              input_nui = lapply(
                seq_along(nui_split_idx), 
                function(i) {
                  method_info[[met_idx]][["nui_unique_nuisance_idx"]][i] %>%
                    # get the nuisance function
                    unique_nuisance[[.]] %>%
                    # get the predicted values of the i-th nuisance function
                    { .[[nui_split_idx[i]]][["predicted.value"]][["test_dataset"]] }
                }) %>%
                # rename according to the nuisance function name
                `names<-`(names(met[["nuisance_fct"]]) %>%
                            `[`( . %in% formalArgs(met[["method"]]) )
                          )
              
              
              c(input_nui, test_dataset) %>%
                # removed unused args
                `[`( names(.) %in% formalArgs(met[["method"]]) ) %>%
                { try(do.call(met[["method"]], .)) } %>%
                # return NA if "met$method" fail
                { if (any(class(.) == "try-error")) NA else . } %>%
                { if (any(is.na(.))) rep(NA, nrow(test_dataset[["X"]])) else . }
            }) %>%
            do.call(cbind, .) %>%
            # compute the mean for each method
            apply(1, mean)
        })
    }) %>%
    # concatenate the result and set the names
    { lapply(seq_along(Balancing_Method), 
             function(i) lapply(., `[[`, i) %>% do.call(cbind, .))
    } %>%
    `names<-`(names(Balancing_Method))
}


crossfit = function(X_tr, Y_tr, Z_tr, X_te, Y_te, Z_te, Balancing_Method,
                    default_nbr_iter=100, default_nbr_fail=NULL) {
  
  # set the default value if not specified
  for (i in seq_along(Balancing_Method)) {
    # default for nbr_iter
    if (is.null(Balancing_Method[[i]]$nbr_iter)) {
      Balancing_Method[[i]]$nbr_iter = default_nbr_iter
    }
    
    # default for max_fail
    if (is.null(Balancing_Method[[i]]$max_fail)) {
      Balancing_Method[[i]]$max_fail = ifelse(
        is.null(default_nbr_fail),
        Balancing_Method[[i]]$nbr_iter,
        default_nbr_fail
      )
    }
    
    # default for nbr_split
    if (is.null(Balancing_Method[[i]]$nbr_split)) {
      Balancing_Method[[i]]$nuisance_fct %>% length
    }
  }
  
  # sort method by "nbr_split"
  Sorted_Balancing_Method = Balancing_Method %>%
    split(x = .,
          f = sapply(., `[[`, "nbr_split")
    )
    
  # cross-fit the methods nbr_iter time
  lapply(
    Sorted_Balancing_Method,
    function(Bal_Met)
    {
      Bal_Met_plus = add_method_info(Bal_Met)
      max_fail = sapply(Bal_Met, `[[`, "max_fail")
      nbr_iter = sapply(Bal_Met, `[[`, "nbr_iter")
      nbr_split = Bal_Met[[1]]$nbr_split
      
      # crossfit the methods to get the right format for output
      output = simple_crossfit(X_tr = X_tr, Y_tr = Y_tr, Z_tr = Z_tr,
                               X_te = X_te, Y_te = Y_te, Z_te = Z_te,
                               nbr_split = nbr_split,
                               nbr_iter = min(nbr_iter),
                               Balancing_Method = Bal_Met_plus)
      
      # set the desired length for each element of output 
      for (i in seq_along(output)) {
        output[[i]] %<>% cbind(matrix(NA, nrow(X_te), nbr_iter[i]-min(nbr_iter)))
      }
      
      # index for which the methods failed
      idx_na = lapply(output, .idx_na_col)
      
      # number of NA per method
      nbr_na = sapply(idx_na, length)
      
      # number of failed exec per method
      nbr_fail = sapply(idx_na, length) - (nbr_iter - min(nbr_iter))
      
      # index of the method to compute
      idx_met = seq_along(Bal_Met)[(nbr_fail < max_fail) & (nbr_na != 0)]
      
      while (length(idx_met) > 0)
      {
        # compute method on new split if they failed
        out = simple_crossfit(
          nbr_iter = min(nbr_na[idx_met]),
          X_tr = X_tr, Y_tr = Y_tr, Z_tr = Z_tr,
          X_te = X_te, Y_te = Y_te, Z_te = Z_te,
          nbr_split = nbr_split,
          Balancing_Method = subset_Balancing_Method(Bal_Met_plus, idx_met))
        
        # save values in output
        for (i in seq_along(out)) {
          j = idx_met[i]
          output[[j]][ , idx_na[[j]][1:min(nbr_na[idx_met])]] = out[[i]]
        }
        
        # update var value
        idx_na = lapply(output, .idx_na_col)
        
        nbr_na = sapply(idx_na, length)
        
        nbr_fail[idx_met] = lapply(out, .idx_na_col) %>%
          sapply(length) + nbr_fail[idx_met]
        
        idx_met = seq_along(Bal_Met)[(nbr_fail < max_fail) & (nbr_na != 0)]
      }
      # print number of failed attempt per method
      print(data.frame(t(nbr_fail), row.names = "nbr fail"))
      
      # compute the median
      output %>%
        lapply(function(met) apply(met, 1, median)) %>%
        return()
    }) %>%
    # concatenate and rename results
    {
      temp = do.call(c, .)
      names(temp) = lapply(., names) %>%
        do.call(c,.)
      temp
    }
}


#### Crossfit Parameters ####


# ## learner
# 
# s_learn = function(mu_interac) {
#   (mu_interac[["1"]] - mu_interac[["0"]]) %>% matrix
# }
# 
# t_learn = function(mu0, mu1) {
#   (mu1 - mu0) %>%  matrix
# }
# 
# x_learn = function(tau0, tau1, g=0.5) {
#   (g*tau0 + (1-g)*tau1) %>% matrix
# }
# 
# 
# ## nuisance function
# 
# # parametric model
# 
# mu_interac_train = function(X, Y, Z) {
#   #no relevant interactions
#   fit <- lrm(Y ~ SEX + rcs(AGE, c(57, 74, 85)) + rcs(SBP, c(130, 160,200)) + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region + Z,
#              data = cbind(X, Y, Z))
#   return(fit)
# }
# 
# mu_interac_pred = function(X, Y, Z, trained.model) {
#   mu = list("0" = NA, "1" = NA)
#   Z = 0
#   mu[["0"]] = predict(trained.model, newdata=cbind(X,Y,Z), type="fitted")
#   Z = 1
#   mu[["1"]] = predict(trained.model, newdata=cbind(X,Y,Z), type="fitted")
#   return(mu)
# }
# 
# mu0_train = function(X, Y, Z) {
#   fit0 <- lrm(Y ~ SEX + rcs(AGE, c(57, 74, 85)) + rcs(SBP, c(130, 160,200)) + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#               data = X,
#               subset = Z==0)
#   return(fit0)
# }
# 
# mu1_train = function(X, Y, Z) {
#   fit1 <- lrm(Y ~ SEX + rcs(AGE, c(57, 74, 85)) + rcs(SBP, c(130, 160,200)) + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#               data = X,
#               subset = Z==1)
#   return(fit1)
# }
# 
# mu0_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type = "fitted")
# }
# 
# mu1_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type = "fitted")
# }
# 
# tau0_train = function(X, Y, Z, mu1) {
#   D0 = mu1[Z==0] - Y[Z==0]
#   tau0 = lm(D0 ~ SEX + rcs(AGE, c(57, 74, 85)) + rcs(SBP, c(130, 160,200)) + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#             data = X[Z==0, ])
#   return(tau0)
# }
# 
# tau1_train = function(X, Y, Z, mu0) {
#   D1 = Y[Z==1] - mu0[Z==1]
#   tau1 = lm(D1 ~ SEX + rcs(AGE, c(57, 74, 85)) + rcs(SBP, c(130, 160,200)) + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#             data = X[Z==1, ])
#   return(tau1)
# }
# 
# tau0_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type="response")
# }
# 
# tau1_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type="response")
# }
# 
# # non parametric model
# 
# mu_interac_rf_train = function(X, Y, Z) {
#   #no relevant interactions
#   fit.rf <- randomForest(Y ~ (SEX + AGE + SBP + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region) * Z,
#                          data = cbind(X, Y, Z),
#                          ntree = 100)
#   return(fit.rf)
# }
# 
# mu_interac_rf_pred = function(X, Y, Z, trained.model) {
#   mu = list("0" = NA, "1" = NA)
#   mu[["0"]] = predict(trained.model, newdata=cbind(X,Y,Z=0), type="response")
#   mu[["1"]] = predict(trained.model, newdata=cbind(X,Y,Z=1), type="response")
#   return(mu)
# }
# 
# 
# 
# mu0_rf_train = function(X, Y, Z) {
#   fit0 <- randomForest(Y ~ SEX + AGE + SBP + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#                        data = X,
#                        subset = Z == 0,
#                        ntree = 100)
#   return(fit0)
# }
# 
# mu1_rf_train = function(X, Y, Z) {
#   fit1 <- randomForest(Y ~ SEX + AGE + SBP + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#                        data = X,
#                        subset = Z == 1,
#                        ntree = 100)
#   return(fit1)
# }
# 
# mu0_rf_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type = "response")
# }
# 
# mu1_rf_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type = "response")
# }
# 
# 
# 
# tau0_rf_train = function(X, Y, Z, mu1) {
#   D0 = mu1 - Y
#   tau0 = randomForest(D0 ~ SEX + AGE + SBP + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#                       data = X,
#                       subset = Z == 0,
#                       ntree = 100)
#   return(tau0)
# }
# 
# tau1_rf_train = function(X, Y, Z, mu0) {
#   D1 = Y - mu0
#   tau1 = randomForest(D1 ~ SEX + AGE + SBP + RDELAY + CTbeforeHosp + CTInfarct + AtrialFib + Asp3d + FaceDef + ArmHandDef + LegFootDef + Dysphasia + Hemianopia + VSDisorder + CerebSigns + OtherDeficit + Conscious + StrokeType + region,
#                       data = X,
#                       subset = Z == 1,
#                       ntree = 100)
#   return(tau1)
# }
# 
# tau0_rf_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type="response")
# }
# 
# tau1_rf_pred = function(X, Y, Z, trained.model) {
#   predict(trained.model, newdata = X, type="response")
# }
# 
# 
# 
# Balancing_Method = list(
#   "S learner" = list(
#     
#     'method' = s_learn,
#     
#     'nuisance_fct' = list(
#       
#       'mu_interac' = list(
#         'train' = mu_interac_train,
#         'predict' = mu_interac_pred)
#     ),
#     
#     'crossfit' = TRUE,
#     'nbr_iter' = 10,
#     'nbr_split' = 5,
#     'max_fail' = 100
#   ),
#   
#   "T learner" = list(
#     
#     'method' = t_learn,
#     
#     'nuisance_fct' = list(
#       
#       'mu0' = list(
#         'train' = mu0_train,
#         'predict' = mu0_pred),
#       
#       'mu1' = list(
#         'train' = mu1_train,
#         'predict' = mu1_pred)
#     ),
#     
#     'crossfit' = TRUE,
#     'nbr_iter' = 10,
#     'nbr_split' = 5,
#     'max_fail' = 100
#   ),
#   
#   "X learner" = list(
#     
#     'method' = x_learn,
#     
#     'nuisance_fct' = list(
#       
#       'mu0' = list(
#         'train' = mu0_train,
#         'predict' = mu0_pred),
#       
#       'mu1' = list(
#         'train' = mu1_train,
#         'predict' = mu1_pred),
#       
#       'tau0' = list(
#         'train' = tau0_train,
#         'predict' = tau0_pred),
#       
#       'tau1' = list(
#         'train' = tau1_train,
#         'predict' = tau1_pred)
#     ),
#     
#     'crossfit' = TRUE,
#     'nbr_iter' = 10,
#     'nbr_split' = 5,
#     'max_fail' = 100
#   )
#   ,
# 
#   "S learner rf" = list(
# 
#     'method' = s_learn,
# 
#     'nuisance_fct' = list(
# 
#       'mu_interac' = list(
#         'train' = mu_interac_rf_train,
#         'predict' = mu_interac_rf_pred)
#     ),
# 
#     'crossfit' = TRUE,
#     'nbr_iter' = 5,
#     'nbr_split' = 5,
#     'max_fail' = 1000
#   ),
# 
#   "T learner rf" = list(
# 
#     'method' = t_learn,
# 
#     'nuisance_fct' = list(
# 
#       'mu0' = list(
#         'train' = mu0_rf_train,
#         'predict' = mu0_rf_pred),
# 
#       'mu1' = list(
#         'train' = mu1_rf_train,
#         'predict' = mu1_rf_pred)
#     ),
# 
#     'crossfit' = TRUE,
#     'nbr_iter' = 5,
#     'nbr_split' = 5,
#     'max_fail' = 1000
#   ),
# 
#   "X learner rf" = list(
# 
#     'method' = x_learn,
# 
#     'nuisance_fct' = list(
# 
#       'mu0' = list(
#         'train' = mu0_rf_train,
#         'predict' = mu0_rf_pred),
# 
#       'mu1' = list(
#         'train' = mu1_rf_train,
#         'predict' = mu1_rf_pred),
# 
#       'tau0' = list(
#         'train' = tau0_rf_train,
#         'predict' = tau0_rf_pred),
# 
#       'tau1' = list(
#         'train' = tau1_rf_train,
#         'predict' = tau1_rf_pred)
#     ),
# 
#     'crossfit' = TRUE,
#     'nbr_iter' = 5,
#     'nbr_split' = 5,
#     'max_fail' = 1000
#   )
# )
#### Crossfit Call & Plot ####


# set.seed(1234)
# system.time(
#   output <- crossfit(X_tr = select(train,-c("deathdep", "aspirin", "group")),
#                      Y_tr = train[ ,"deathdep"],
#                      Z_tr = train[ ,"aspirin"],
#                      X_te = select(test,-c("deathdep", "aspirin", "group")),
#                      Y_te = test[ ,"deathdep"],
#                      Z_te = test[ ,"aspirin"],
#                      Balancing_Method)
# )

#### Plots ####

### heatmap

heatmap = function(data, itr1, itr2, covariates, target = "MCC") {
  if (!(target %in% c("MCC", "Kappa")))
    stop("'target' must be either 'MCC' or 'Kappa'")
  data = cbind(data[covariates],
               "n11" =  itr1 &  itr2,
               "n10" =  itr1 & !itr2,
               "n01" = !itr1 &  itr2,
               "n00" = !itr1 & !itr2,
               "n"   =  1)
  lapply(
    covariates,
    function(c1)
    {
      lapply(
        covariates,
        function(c2)
        {
          df = aggregate(
            x = data[c("n11", "n10", "n01", "n00", "n")],
            by = data[c(c1, c2)],
            FUN = sum)
          
          # compute MCC / Phi coefficient
          if (target == "MCC")
            df %<>% mutate(
              "MCC" = (n11n00 - n10n01) / exp(
                0.5 * (log(n11+n10) + log(n01+n00) + log(n10+n00) + log(n11+n01))
              ))
          
          # compute Cohen's kappa coefficient
          if (target == "Kappa")
            df %<>% mutate(
              "Kappa" = 1 - n * (
                (n - n11 - n00) / (n^2 - (n11+n10)(n11+n01) - (n01+n00)(n10+n00))
              ))
          
          df %<>% [(-(3:7)) # remove variables 'n11', 'n10', 'n01', 'n00', 'n'
                   df[ ,1] = paste(c1,"=",df[ ,1])
                   df[ ,2] = paste(c2,"=",df[ ,2])
                   names(df)[1:2] = c("cov1", "cov2")
                   return(df)
        })
    }) %>%
    do.call("c", .) %>%
    do.call("rbind", .) %>%
    ggplot( aes(x = cov1, y = cov2, fill = .data[[target]]) ) +
    geom_tile() +
    labs(fill = target) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

# covariates = select(test, -c("AGE", "RDELAY", "SBP", "group")) %>% names
# heatmap(test, output[["X learner"]]<0, output[["T learner"]]<0, covariates)
###Scatter plot

ite_scatter_plot = function(data=NULL, ite1, ite2, label=c("ite1", "ite2"),
                            title = "ITE Distribution") {
  if (!is.null(data) & typeof(data) != "list" )
    stop("'data' must be a list")
  if (is.null(data)) data = list(ite1, ite2) %>% `names<-`(label)
  if (any(sapply(data, is.character)))
    stop("if 'data' is not given then 'ite1' and 'ite2' must be numeric")
  if (is.character(c(ite1,ite2))) label = c(ite1, ite2)
  
  ite1 = data[[label[1]]]
  ite2 = data[[label[2]]]
  
  window_size = sapply(data, range) %>%
    abs %>% max %>% c(-., .) * 1.05
  
  ggplot() +
    xlim(window_size) +
    ylim(window_size) +
    geom_point(aes(x=ite1, y=ite2), size=1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(label[1]) +
    ylab(label[2]) +
    ggtitle(label = title)
}


# ite_scatter_plot(output, output[["T learner"]], output[["X learner"]], c("T learner", "X learner"))

# density plot

# ggplot() +
#   geom_density(aes(x=output[["S learner"]], color="S learner")) +
#   geom_density(aes(x=output[["T learner"]], color="T learner")) +
#   geom_density(aes(x=output[["X learner"]], color="X learner")) +
#   geom_density(aes(x=output[["S learner rf"]], color="S learner rf")) +
#   geom_density(aes(x=output[["T learner rf"]], color="T learner rf")) +
#   geom_density(aes(x=output[["X learner rf"]], color="X learner rf")) +
#   geom_vline(xintercept = 0) +
#   xlab("ITE") +
#   ylab("Density") +
#   ggtitle("ITE Density")




